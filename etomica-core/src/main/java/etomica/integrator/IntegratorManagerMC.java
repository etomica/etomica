/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.mcmove.*;
import etomica.util.Arrays;
import etomica.util.IEvent;
import etomica.util.IEventManager;
import etomica.util.random.IRandom;

/**
 * Integrator manages other Integrators which either act on a Box, or manage
 * still other Integrators.  Each step, this class either performs global MC 
 * moves between the Integrators or runs the sub-integrators.
 * 
 * @author David Kofke and Andrew Schultz
 */

public class IntegratorManagerMC extends Integrator {

    private static final long serialVersionUID = 2L;
    protected final IEventManager eventManager;
    protected final IRandom random;
    private final IEvent trialEvent;
    private final IEvent acceptedEvent, rejectedEvent;
    protected double globalMoveProbability;
    protected MCMoveManager moveManager;
    protected Integrator[] integrators;
    protected int nIntegrators;
    protected double temperature;
    private double globalMoveInterval;
    
    public IntegratorManagerMC(IRandom random) {
        super();
        this.random = random;
        integrators = new Integrator[0];
        setGlobalMoveInterval(2);
        moveManager = new MCMoveManager(random);
        eventManager = new MCMoveEventManager();
        trialEvent = new MCMoveTrialInitiatedEvent(moveManager);
        acceptedEvent = new MCMoveTrialCompletedEvent(moveManager, true);
        rejectedEvent = new MCMoveTrialCompletedEvent(moveManager, false);
    }

    protected void setup() {
        ConfigurationOverlapException overlapException = null;
        for(int i=0; i<nIntegrators; i++) {
            try {
                integrators[i].reset();
            }
            catch (ConfigurationOverlapException e) {
                if (overlapException == null) {
                    overlapException = e;
                }
            }
        }
        super.setup();
        if (overlapException != null) {
            throw overlapException;
        }
    }
    
    /**
     * Causes recalculation of move frequencies and zero of selection counts for
     * moves.  Resets this integrator and passes on the reset to all managed integrators.
     */
    public void reset() throws ConfigurationOverlapException {
        super.reset();

        moveManager.recomputeMoveFrequencies();
        ConfigurationOverlapException overlapException = null;
        for(int i=0; i<integrators.length; i++) {
            try {
                integrators[i].reset();
            }
            catch (ConfigurationOverlapException e) {
                if (overlapException == null) {
                    overlapException = e;
                }
            }
        }
        if (overlapException != null) {
            throw overlapException;
        }
    }

    public void addIntegrator(Integrator integrator){
        integrators = (Integrator[])Arrays.addObject(integrators,integrator);
        nIntegrators++;
    }

    /**
     * Removes the given integrator from the list of integrators.  Returns
     * false if the given integrator was not handled by this integrator.
     */
    public boolean removeIntegrator(Integrator integrator) {
        integrators = (Integrator[])Arrays.removeObject(integrators,integrator);
        if (nIntegrators == integrators.length) {
            return false;
        }
        nIntegrators--;
        return true;
    }

    public Integrator[] getIntegrators() {
        return integrators.clone();
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
     * Performs a Monte Carlo trial that attempts to swap the configurations
     * between two "adjacent" boxes, or instructs all integrators to perform
     * a single doStep.
     */
    protected void doStepInternal() {
        if(random.nextDouble() < globalMoveProbability) {
            doGlobalMoves();
        } else {
            for(int i=0; i<nIntegrators; i++) {
                integrators[i].doStep();
            }
        }
    }

    /**
     * Method to select and perform an elementary Monte Carlo move. The type of
     * move performed is chosen from all MCMoves that have been added to the
     * integrator. Each MCMove has associated with it a (unnormalized)
     * frequency, which when weighed against the frequencies given the other
     * MCMoves, determines the likelihood that the move is selected. After
     * completing move, fires an MCMove event if there are any listeners.
     */
    protected void doGlobalMoves() {
        //select the move
        MCMove move = moveManager.selectMove();
        if (move == null)
            return;

        //perform the trial
        //returns false if the trial cannot be attempted; for example an
        // atom-displacement trial in a box with no molecules
        if (!move.doTrial())
            return;

        //notify any listeners that move has been attempted
        eventManager.fireEvent(trialEvent);

        //decide acceptance
        double chi = move.getChi(temperature);
        if (chi == 0.0 || (chi < 1.0 && chi < random.nextDouble())) {
            //reject
            move.rejectNotify();
            //notify listeners of outcome
            eventManager.fireEvent(rejectedEvent);
            move.getTracker().updateCounts(false, chi);
        } else {
            //accept
            move.acceptNotify();
            //notify listeners of outcome
            eventManager.fireEvent(acceptedEvent);
            move.getTracker().updateCounts(true, chi);
        }
    }

    public IEventManager getMoveEventManager() {
        return eventManager;
    }

    /**
     * Accessor method for the average interval between box-swap trials.
     */
    public double getGlobalMoveInterval() {
        return globalMoveInterval;
    }

    /**
     * Sets the average interval between box-swap trials.  With each
     * call to doStep of this integrator, there will be a probability of
     * 1/globalMoveInterval that a swap trial will be attempted.  A swap is attempted
     * for only one pair of boxs.  Default is 2.
     */
    public void setGlobalMoveInterval(double newGlobalMoveInterval) {
        if(newGlobalMoveInterval <= 0) {
            throw new IllegalArgumentException("global move interval must be positive");
        }
        globalMoveInterval = newGlobalMoveInterval;
        if (Double.isInfinite(globalMoveInterval)) {
            globalMoveProbability = 0;
        }
        else {
            globalMoveProbability = 1.0/globalMoveInterval;
        }
    }

    public double getTemperature() {
        return temperature;
    }

    /**
     * Sets the temperature for this integrator if that's relevant
     * @param temperature
     */
    public void setTemperature(double temperature) {
        this.temperature = temperature;
    }
}
