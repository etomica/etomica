/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.api.IRandom;
import etomica.simulation.Simulation;
import etomica.integrator.mcmove.MCMove;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.integrator.mcmove.MCMoveEventManager;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.integrator.mcmove.MCMoveTrialCompletedEvent;
import etomica.integrator.mcmove.MCMoveTrialFailedEvent;
import etomica.integrator.mcmove.MCMoveTrialInitiatedEvent;
import etomica.util.IEvent;
import etomica.util.IEventManager;

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

public class IntegratorMC extends IntegratorBox {

    public IntegratorMC(Simulation sim, PotentialMaster potentialMaster) {
        this(potentialMaster, sim.getRandom(), 1.0);
    }
    
	/**
	 * Constructs integrator and establishes PotentialMaster instance that
	 * will be used by moves to calculate the energy.
	 */
	public IntegratorMC(PotentialMaster potentialMaster, IRandom random, double temperature) {
		super(potentialMaster,temperature);
        this.random = random;
		setIsothermal(true); //has no practical effect, but sets value of
		// isothermal to be consistent with way integrator
		// is sampling
        moveManager = new MCMoveManager(random);
        moveEventManager = new MCMoveEventManager();
        trialEvent = new MCMoveTrialInitiatedEvent(moveManager);
        trialFailedEvent = new MCMoveTrialFailedEvent(moveManager);
        acceptedEvent = new MCMoveTrialCompletedEvent(moveManager, true);
        rejectedEvent = new MCMoveTrialCompletedEvent(moveManager, false);
	}

    /**
     * @return Returns the moveManager.
     */
    public MCMoveManager getMoveManager() {
        return moveManager;
    }

    /**
     * @param newMoveManager The moveManager to set.
     */
    public void setMoveManager(MCMoveManager newMoveManager) {
        moveManager = newMoveManager;
    }
    
    /**
     * Invokes superclass method and informs all MCMoves about the new box.
     * Moves are not notified if they have a number of boxs different from
     * the number of boxs handled by the integrator.
     */
    public void setBox(Box p) {
    	super.setBox(p);
    	moveManager.setBox(p);
    }

    /**
     * Method to select and perform an elementary Monte Carlo move. The type of
     * move performed is chosen from all MCMoves that have been added to the
     * integrator. Each MCMove has associated with it a (unnormalized)
     * frequency, which when weighed against the frequencies given the other
     * MCMoves, determines the likelihood that the move is selected. After
     * completing move, fires an MCMove event if there are any listeners.
     */
    public void doStepInternal() {
    	//select the move
    	MCMoveBox move = (MCMoveBox)moveManager.selectMove();
    	if (move == null)
    		return;

    	//perform the trial
    	//returns false if the trial cannot be attempted; for example an
    	// atom-displacement trial in a box with no molecules
    	if (!move.doTrial()) {
    	    moveEventManager.fireEvent(trialFailedEvent);
    		return;
    	}

        //notify any listeners that move has been attempted
		moveEventManager.fireEvent(trialEvent);

    	//decide acceptance
    	double chi = move.getA() * Math.exp(move.getB()/temperature);
    	if (chi == 0.0 || (chi < 1.0 && chi < random.nextDouble())) {//reject
            if (dodebug) {
                System.out.println(stepCount+" move "+move+" rejected");
            }
            move.getTracker().updateCounts(false, chi);
    		move.rejectNotify();
            //notify listeners of outcome
            moveEventManager.fireEvent(rejectedEvent);
    	} else {
            if (dodebug) {
                System.out.println(stepCount+" move "+move+" accepted");
            }
            move.getTracker().updateCounts(true, chi);
    		move.acceptNotify();
    		currentPotentialEnergy += move.energyChange();
            //notify listeners of outcome
            moveEventManager.fireEvent(acceptedEvent);
    	}
    }
    
    public void notifyEnergyChange(double energyChange) {
        currentPotentialEnergy += energyChange;
    }

    /**
     * Causes recalculation of move frequencies and zero of selection counts for
     * moves.
     */
    public void reset() {
        super.reset();
        moveManager.recomputeMoveFrequencies();
    }

    /**
     * Adds a listener that will be notified when a MCMove trial is attempted
     * and when it is completed.
     */
    public IEventManager getMoveEventManager() {
        return moveEventManager;
    }

    private static final long serialVersionUID = 2L;
    protected final IRandom random;
    protected MCMoveManager moveManager;
    protected final IEventManager moveEventManager;
    private final IEvent trialEvent, trialFailedEvent;
    private final IEvent acceptedEvent, rejectedEvent;
    public static boolean dodebug;
}
