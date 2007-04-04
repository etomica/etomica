package etomica.integrator;

import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.mcmove.MCMove;
import etomica.integrator.mcmove.MCMoveEventManager;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.integrator.mcmove.MCMoveTrialCompletedEvent;
import etomica.integrator.mcmove.MCMoveTrialInitiatedEvent;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.util.Arrays;
import etomica.util.IRandom;

/**
 * Integrator manages other Integrators which either act on a Phase, or manager 
 * still other Integrators.  Each step, this class either performs global MC 
 * moves between the Integrators or runs the sub-integratos.
 * 
 * @author David Kofke and Andrew Schultz
 */

public class IntegratorManagerMC extends Integrator {

    public IntegratorManagerMC(Simulation sim) {
        this(sim.getPotentialMaster(), sim.getRandom());
    }
    
    public IntegratorManagerMC(PotentialMaster potentialMaster, IRandom random) {
        super(potentialMaster);
        this.random = random;
        integrators = new Integrator[0];
        intervalEvents = new IntegratorIntervalEvent[0];
        setGlobalMoveInterval(2);
        moveManager = new MCMoveManager(random);
        eventManager = new MCMoveEventManager();
        trialEvent = new MCMoveTrialInitiatedEvent(moveManager);
        acceptedEvent = new MCMoveTrialCompletedEvent(moveManager, true);
        rejectedEvent = new MCMoveTrialCompletedEvent(moveManager, false);
    }

    protected void setup() throws ConfigurationOverlapException {
        ConfigurationOverlapException overlapException = null;
        for(int i=0; i<nIntegrators; i++) {
            try {
                integrators[i].initialize();
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
        
        super.reset();
    }

    public void addIntegrator(Integrator integrator){
        integrators = (Integrator[])Arrays.addObject(integrators,integrator);
        intervalEvents = (IntegratorIntervalEvent[])Arrays.addObject(intervalEvents,new IntegratorIntervalEvent(integrator, 1));
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
        return (Integrator[])integrators.clone();
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
     * Fires non-interval event for this integrator, then instructs
     * each sub-integrator to fire event.
     */
    public void fireNonintervalEvent(IntegratorNonintervalEvent ie) {
        super.fireNonintervalEvent(ie);
        if (ie.type() == IntegratorNonintervalEvent.DONE) {
            for(int i=0; i<nIntegrators; i++) {
                integrators[i].fireNonintervalEvent(new IntegratorNonintervalEvent(integrators[i], IntegratorNonintervalEvent.DONE));
            }
        }
    }
    
    /**
     * Performs a Monte Carlo trial that attempts to swap the configurations
     * between two "adjacent" phases, or instructs all integrators to perform
     * a single doStep.
     */
    public void doStep() {
        if(random.nextDouble() < globalMoveProbability) {
            doGlobalMoves();
        } else {
            for(int i=0; i<nIntegrators; i++) {
                integrators[i].doStep();
                integrators[i].fireIntervalEvent(intervalEvents[i]);
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
        // atom-displacement trial in a phase with no molecules
        if (!move.doTrial())
            return;

        //notify any listeners that move has been attempted
        eventManager.fireEvent(trialEvent);

        //decide acceptance
        double chi = move.getA();
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

    public MCMoveEventManager getMoveEventManager() {
        return eventManager;
    }
    
    /**
     * Sets the average interval between phase-swap trials.  With each 
     * call to doStep of this integrator, there will be a probability of
     * 1/globalMoveInterval that a swap trial will be attempted.  A swap is attempted
     * for only one pair of phases.  Default is 2.
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
    
    /**
     * Accessor method for the average interval between phase-swap trials.
     */
    public double getGlobalMoveInterval() {return globalMoveInterval;}
    
    private static final long serialVersionUID = 2L;
    private double globalMoveInterval;
    protected double globalMoveProbability;
    protected MCMoveManager moveManager;
    protected final MCMoveEventManager eventManager;
    protected Integrator[] integrators;
    protected IntegratorIntervalEvent[] intervalEvents;
    protected int nIntegrators;
    private final MCMoveTrialInitiatedEvent trialEvent;
    private final MCMoveTrialCompletedEvent acceptedEvent, rejectedEvent;
    protected final IRandom random;
}
