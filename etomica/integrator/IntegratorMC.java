package etomica.integrator;

import etomica.Atom;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Integrator;
import etomica.Phase;
import etomica.PotentialMaster;
import etomica.Simulation;
import etomica.SimulationEventManager;
import etomica.integrator.mcmove.MCMoveEvent;
import etomica.integrator.mcmove.MCMoveListener;

/**
 * Integrator to perform Metropolis Monte Carlo sampling. Works with a set of
 * MCMove instances that are added to the integrator. A step performed by the
 * integrator consists of selecting a MCMove from the set, performing the trial
 * defined by the MCMove, and deciding acceptance of the trial using information
 * from the MCMove.
 * 
 * @author David Kofke
 * @see MCMove
 */

public class IntegratorMC extends Integrator implements EtomicaElement {

	/**
	 * Constructs integrator and establishes PotentialMaster instance that
	 * will be used by moves to calculate the energy.
	 */
	public IntegratorMC(PotentialMaster potentialMaster) {
		super(potentialMaster);
		setIsothermal(true); //has no practical effect, but sets value of
		// isothermal to be consistent with way integrator
		// is sampling
	}

	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo("General Monte Carlo simulation");
		return info;
	}

	/**
	 * Sets moves in given array to be integrator's set of moves, deleting any
	 * existing moves.
	 */
	public void setMCMoves(MCMove[] moves) {
		firstMoveLink = null;
		moveCount = 0;
		for (int i = 0; i < moves.length; i++) {
			addMCMove(moves[i]);
		}
	}

	/**
	 * Constructs and returns array of all moves added to the integrator.
	 */
	public MCMove[] getMCMoves() {
		MCMove[] moves = new MCMove[moveCount];
		int i = 0;
		for (MCMoveLinker link = firstMoveLink; link != null; link = link.nextLink) {
			moves[i++] = link.move;
		}
		return moves;
	}

	/**
	 * Adds the given MCMove to the set of moves performed by the integrator and
	 * recalculates move frequencies.
	 */
	public void addMCMove(MCMove move) {
		//make sure move wasn't added already
		for (MCMoveLinker link = firstMoveLink; link != null; link = link.nextLink) {
			if (move == link.move)
				return;
		}
		if (firstMoveLink == null) {
			firstMoveLink = new MCMoveLinker(move);
            lastMoveLink = firstMoveLink;
		} else {
			lastMoveLink.nextLink = new MCMoveLinker(move);
			lastMoveLink = lastMoveLink.nextLink;
		}
		if(move.phases.length == phase.length) move.setPhase(phase);
		move.setTemperature(temperature);
		moveCount++;
		recomputeMoveFrequencies();
	}

	/**
	 * Invokes superclass method and informs all MCMoves about the new phase.
     * Moves are not notified if they have a number of phases different from
     * the number of phases handled by the integrator.
	 */
	public boolean addPhase(Phase p) {
		if (!super.addPhase(p))
			return false;
		for (MCMoveLinker link = firstMoveLink; link != null; link = link.nextLink) {
			if(link.move.phases.length == phase.length) link.move.setPhase(phase);
		}
		return true;
	}

	/**
	 * Selects a MCMove instance from among those added to the integrator, with
	 * probability in proportion to the frequency value assigned to the move.
	 */
	protected MCMove selectMove() {
		if (firstMoveLink == null)
			return null;
		int i = Simulation.random.nextInt(frequencyTotal);
		MCMoveLinker link = firstMoveLink;
		while ((i -= link.fullFrequency) >= 0) {
			link = link.nextLink;
		}
		link.selectionCount++;
		return link.move;
	}

	/**
	 * Method to select and perform an elementary Monte Carlo move. The type of
	 * move performed is chosen from all MCMoves that have been added to the
	 * integrator. Each MCMove has associated with it a (unnormalized)
	 * frequency, which when weighed against the frequencies given the other
	 * MCMoves, determines the likelihood that the move is selected. After
	 * completing move, fires an MCMove event if there are any listeners.
	 */
	public void doStep() {
		//select the move
		MCMove move = selectMove();
		if (move == null)
			return;

		//perform the trial
		//returns false if the trial cannot be attempted; for example an
		// atom-displacement trial in a phase with no molecules
		if (!move.doTrial())
			return;

		//notify any listeners that move has been attempted
		if (eventManager != null) { //consider using a final boolean flag that
			// is set in constructor
			event.mcMove = move;
			event.isTrialNotify = true;
			eventManager.fireEvent(event);
		}

		//decide acceptance
		double lnChi = move.lnTrialRatio() + move.lnProbabilityRatio();
        double chi = lnChi == -Double.POSITIVE_INFINITY ? 0.0 : 
                                (lnChi > 0.0 ? 1.0 : Math.exp(lnChi));
		if (chi == 0.0 || (chi < 1.0 && chi < Simulation.random.nextDouble())) {//reject
			move.rejectNotify();
			event.wasAccepted = false;
		} else {
			move.acceptNotify();
			event.wasAccepted = true;
            for (int i=0; i<phase.length; i++) {
                currentPotentialEnergy[i] += move.energyChange(phase[i]);
            }
		}

		//notify listeners of outcome
		if (eventManager != null) {
			event.isTrialNotify = false;
			eventManager.fireEvent(event);
		}

		move.updateCounts(event.wasAccepted, chi, equilibrating);
	}

	/**
	 * Sets the partial, unnormalized frequency for performing the given move,
	 * relative to the other moves that have been added to the integrator. If
	 * the perParticleFrequency flag for the move is true, the full frequency is
	 * determined by multiplying this partial frequency by the number of
	 * molecules in the phases affected by the integrator when this method (or
	 * setPerParticleFrequency) is invoked; otherwise the full frequency equals
	 * this partial frequency. <br>
	 * Each move is performed (on average) an amount in proportion to the full
	 * frequency. Moves having the same full frequency are performed with equal
	 * likelihood. <br>
	 * Default value of the (partial) frequency is 100, which is (for example)
	 * the nominal value for MCMoveAtom. <br>
	 */

	public void setFrequency(MCMove move, int frequency) {
		MCMoveLinker link = null;
		for (link = firstMoveLink; link != null; link = link.nextLink) {
			if (link.move == move)
				break;
		}
		if (link == null)
			throw new RuntimeException(
					"Attempt to change frequency of MCMove that is not managed by integrator");
		link.frequency = frequency;
		recomputeMoveFrequencies();
	}

	/**
	 * Indicates if frequency for the given move indicates full frequency, or
	 * frequency per particle. If per particle, frequency assigned to move is
	 * multiplied by the number of particles affected by the move at time this
	 * method (or setFrequency) is called (full frequency is not otherwise
	 * regularly updated for changing particle numbers). <br>
	 * 
	 * @see #setFrequency
	 */
	public final void setPerParticleFrequency(MCMove move,
			boolean isPerParticleFrequency) {
		MCMoveLinker link = null;
		for (link = firstMoveLink; link != null; link = link.nextLink) {
			if (link.move == move)
				break;
		}
		if (link == null)
			throw new RuntimeException(
					"Attempt to change frequency of MCMove that is not managed by integrator");
		link.perParticleFrequency = isPerParticleFrequency;
		recomputeMoveFrequencies();
	}

	/**
	 * Recomputes all the move frequencies.
	 */
	protected void recomputeMoveFrequencies() {
		frequencyTotal = 0;
		for (MCMoveLinker link = firstMoveLink; link != null; link = link.nextLink) {
			link.resetFullFrequency();
			frequencyTotal += link.fullFrequency;
		}
	}

	/**
	 * Sets the temperature for this integrator and all the MCMove instances it
	 * currently holds.
	 */
	public void setTemperature(double temperature) {
		super.setTemperature(temperature);
		for (MCMoveLinker link = firstMoveLink; link != null; link = link.nextLink) {
			link.move.setTemperature(temperature);
		}
	}

    /**
     * Causes recalculation of move frequencies and zero of selection counts for
     * moves.
     */
    public void reset() {
        recomputeMoveFrequencies();
        for (MCMoveLinker link = firstMoveLink; link != null; link = link.nextLink) {
            link.selectionCount = 0;
        }
        super.reset();
    }

	/**
	 * Adds a listener that will be notified when a MCMove trial is attempted
	 * and when it is completed.
	 */
	public void addMCMoveListener(MCMoveListener listener) {
		if (eventManager == null)
			eventManager = new SimulationEventManager();
		eventManager.addListener(listener);
	}

	/**
	 * Removes the given listener.
	 */
	public void removeMCMoveListener(MCMoveListener listener) {
		if (eventManager == null)
			return; //should define an exception
		eventManager.removeListener(listener);
		if (eventManager.listenerCount() == 0)
			eventManager = null;
	}

	/**
	 * Returns null.
	 */
	public Object makeAgent(Atom a) {
		return null;
	}

	private MCMoveLinker firstMoveLink, lastMoveLink;
	private int frequencyTotal;
	private int moveCount;
	protected SimulationEventManager eventManager;
	protected final MCMoveEvent event = new MCMoveEvent(this);

	/**
	 * Linker used to construct linked-list of MCMove instances
	 */
	private static class MCMoveLinker {
		int frequency, fullFrequency;
		final MCMove move;
		boolean perParticleFrequency;
		int selectionCount;
		MCMoveLinker nextLink;

		MCMoveLinker(MCMove move) {
			this.move = move;
			frequency = move.getNominalFrequency();
			perParticleFrequency = move.isNominallyPerParticleFrequency();
		}

		/**
		 * Updates the full frequency based on the current value of the
		 * frequency, the status of the perParticleFrequency flag, and the
		 * current number of molecules in the phases affected by the move.
		 */
		void resetFullFrequency() {
			fullFrequency = frequency;
            Phase[] phases = move.getPhase();
			if (perParticleFrequency && phases != null) {
				int mCount = 0;
				for (int i = 0; i < phases.length; i++)
					if (phases[i] != null)
						mCount += phases[i].moleculeCount();
				fullFrequency *= mCount;
			}
		}

	}
}//end of IntegratorMC
