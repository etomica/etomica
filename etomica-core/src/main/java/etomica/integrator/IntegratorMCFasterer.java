/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.box.Box;
import etomica.integrator.mcmove.*;
import etomica.potential.compute.PotentialCompute;
import etomica.simulation.Simulation;
import etomica.util.EventManager;
import etomica.util.random.IRandom;

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

public class IntegratorMCFasterer extends IntegratorBoxFasterer {

    public static boolean dodebug;
    protected final IRandom random;
    protected final EventManager<MCMoveEvent> moveEventManager;
    private final MCMoveEvent trialEvent, trialFailedEvent;
    private final MCMoveEvent acceptedEvent, rejectedEvent;
    protected MCMoveManager moveManager;

    /**
     * @param sim             Simulation where this integrator is used
     * @param potentialCompute PotentialMaster instance used by moves to calculate the energy
     */

    public IntegratorMCFasterer(Simulation sim, PotentialCompute potentialCompute, Box box) {
        this(potentialCompute, sim.getRandom(), 1.0, box);
    }

    /**
     * @param potentialCompute PotentialMaster instance used by moves to calculate the energy
     * @param random          random number generator used to select moves and decide acceptance
     * @param temperature     temperature of the ensemble
     */
    public IntegratorMCFasterer(PotentialCompute potentialCompute, IRandom random, double temperature, Box box) {
        super(potentialCompute, temperature, box);
        this.random = random;
        setIsothermal(true); //has no practical effect, but sets value of
        // isothermal to be consistent with way integrator
        // is sampling
        moveManager = new MCMoveManager(random);
        moveEventManager = new EventManager<>();
        trialEvent = new MCMoveTrialInitiatedEvent(moveManager);
        trialFailedEvent = new MCMoveTrialFailedEvent(moveManager);
        acceptedEvent = new MCMoveTrialCompletedEvent(moveManager, true);
        rejectedEvent = new MCMoveTrialCompletedEvent(moveManager, false);
        moveManager.setBox(box);
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
     * Method to select and perform an elementary Monte Carlo move and decide acceptance. The type of
     * move performed is chosen from all MCMoves that have been added to the
     * integrator. After completing move, fires an MCMove event if there are any listeners.
     */
    protected void doStepInternal() {
        //select the move
        MCMoveBox move = (MCMoveBox) moveManager.selectMove();
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
        double chi = move.getChi(temperature);
        if (chi == 0.0 || (chi < 1.0 && chi < random.nextDouble())) {//reject
            if (dodebug) {
                System.out.println(stepCount + " move " + move + " rejected " + chi);
            }
            move.getTracker().updateCounts(false, chi);
            move.rejectNotify();
            //notify listeners of outcome
            moveEventManager.fireEvent(rejectedEvent);
        } else {
            if (dodebug) {
                System.out.println(stepCount + " move " + move + " accepted " + chi);
            }
            move.getTracker().updateCounts(true, chi);
            move.acceptNotify();
            currentPotentialEnergy += move.energyChange();
            //notify listeners of outcome
            moveEventManager.fireEvent(acceptedEvent);
        }
    }

    /**
     * Notifies the IntegratorMC that the energy changed allowing
     * the internal potential energy field to be updated.
     *
     * @param energyChange Change in the energy
     */
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
     * @return moveEventManager that fires move events
     */
    public EventManager<MCMoveEvent> getMoveEventManager() {
        return moveEventManager;
    }
}
