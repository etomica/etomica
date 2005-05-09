package etomica.action;

import java.util.Iterator;

import etomica.Integrator;
import etomica.Phase;
import etomica.Simulation;
import etomica.data.DataAccumulator;

/**
 * Action that invokes reset method of all registered simulation elements,
 * effectively initializing the entire simulation.
 */
 
 /* History of changes
  * 7/03/02 (DAK/SKK) Modified to loop through all elements, using reset method (new to SimulationElement).
  */

public final class SimulationRestart extends SimulationActionAdapter {
    
    public SimulationRestart(Simulation sim) {
        super("Reset");
        setSimulation(sim);
    }
        
    /**
     * Resets phases, integrators, and accumulators.
     */
    public void actionPerformed() {

        for(Iterator iter=simulation.getPhaseList().iterator(); iter.hasNext(); ) {
            Phase phase = (Phase)iter.next();
            phase.reset();
        }
        
        for(Iterator iter=simulation.getIntegratorList().iterator(); iter.hasNext(); ) {
            Integrator integrator = (Integrator)iter.next();
            if(integrator.isInitialized()) integrator.initialize();
        }
        
        for(Iterator iter=simulation.getDataAccumulatorList().iterator(); iter.hasNext(); ) {
            DataAccumulator dataAccumulator = (DataAccumulator)iter.next();
            dataAccumulator.reset();
        }
    }
    
}
        