package etomica.integrator;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.mcmove.MCMoveEvent;
import etomica.integrator.mcmove.MCMoveEventManager;
import etomica.integrator.mcmove.MCMoveListener;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;

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

public class IntegratorMC extends IntegratorPhase implements EtomicaElement {

    public IntegratorMC(Simulation sim) {
        this(sim.potentialMaster,sim.getDefaults().temperature);
    }
    
	/**
	 * Constructs integrator and establishes PotentialMaster instance that
	 * will be used by moves to calculate the energy.
	 */
	public IntegratorMC(PotentialMaster potentialMaster, double temperature) {
		super(potentialMaster,temperature);
		setIsothermal(true); //has no practical effect, but sets value of
		// isothermal to be consistent with way integrator
		// is sampling
        moveManager = new MCMoveManager();
	}

	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo("General Monte Carlo simulation");
		return info;
	}

    /**
     * @return Returns the moveManager.
     */
    public MCMoveManager getMoveManager() {
        return moveManager;
    }

    /**
     * @param moveManager The moveManager to set.
     */
    public void setMoveManager(MCMoveManager newMoveManager) {
        moveManager = newMoveManager;
    }

    /**
     * Invokes superclass method and informs all MCMoves about the new phase.
     * Moves are not notified if they have a number of phases different from
     * the number of phases handled by the integrator.
     */
    public void setPhase(Phase p) {
    	super.setPhase(p);
    	moveManager.setPhase(p);
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
    	MCMove move = moveManager.selectMove();
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
    		currentPotentialEnergy += move.energyChange(phase);
    	}

    	//notify listeners of outcome
    	if (eventManager != null) {
    		event.isTrialNotify = false;
    		eventManager.fireEvent(event);
    	}

    	move.updateCounts(event.wasAccepted, chi, equilibrating);
    }

    /**
     * Sets the temperature for this integrator and all the MCMove instances it
     * currently holds.
     */
    public void setTemperature(double temperature) {
    	super.setTemperature(temperature);
    	// moveManager will be null when called from superclass constructor
    	if (moveManager == null) {
    		moveManager.setTemperature(temperature);
    	}
    }

    /**
     * Causes recalculation of move frequencies and zero of selection counts for
     * moves.
     */
    public void reset() throws ConfigurationOverlapException {
        moveManager.recomputeMoveFrequencies();
        super.reset();
    }

    /**
     * Adds a listener that will be notified when a MCMove trial is attempted
     * and when it is completed.
     */
    public void addMCMoveListener(MCMoveListener listener) {
    	if (eventManager == null)
    		eventManager = new MCMoveEventManager();
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

    protected MCMoveManager moveManager;
    protected MCMoveEventManager eventManager;
	protected final MCMoveEvent event = new MCMoveEvent(this);

}
