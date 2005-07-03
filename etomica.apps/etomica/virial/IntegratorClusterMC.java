package etomica.virial;

import etomica.PotentialMaster;
import etomica.Simulation;
import etomica.integrator.IntegratorMC;
import etomica.integrator.MCMove;

/**
 * Integrator appropriate for cluster simulations.  Moves are assumed to
 * be instances of MCMoveCluster.  The sampling weight bias of the current 
 * configuration is tracked so property averages can be unbiased efficiently 
 * (without recalculating the bias).  
 */

/*
 * Created on Sep 10, 2004
 */
public class IntegratorClusterMC extends IntegratorMC {

	protected double weight;
	
    public IntegratorClusterMC(PotentialMaster potentialMaster) {
        super(potentialMaster);
        weight = 1.0;
    }

	/**
     * Method to select and perform an elementary Monte Carlo cluster move.  
     */
    public void doStep() {
        //select the move
        MCMove move = selectMove();
        if(move == null) return;
        
        //perform the trial
        //returns false if the trial cannot be attempted; for example an atom-displacement trial in a phase with no molecules
        if(!move.doTrial()) return;
        
        //notify any listeners that move has been attempted
        if(eventManager != null) { //consider using a final boolean flag that is set in constructor
            event.mcMove = move;
            event.isTrialNotify = true;
            eventManager.fireEvent(event);
        }
        
        //decide acceptance
        double chi = ((MCMoveCluster)move).trialRatio() * ((MCMoveCluster)move).probabilityRatio();
        if(chi == 0.0 || (chi < 1.0 && chi < Simulation.random.nextDouble())) {//reject
            move.rejectNotify();
            event.wasAccepted = false;
        } else {
            move.acceptNotify();
            event.wasAccepted = true;
            weight *= chi;
        }

        //notify listeners of outcome
        if(eventManager != null) { //consider using a final boolean flag that is set in constructor
            event.isTrialNotify = false;
            eventManager.fireEvent(event);
        }
        
        move.updateCounts(event.wasAccepted,chi,isEquilibrating());
    }
    
    public double getWeight() {
    	return weight;
    }
    
    public void reset() {
    	super.reset();
    	weight = 1.0;
    }
}
