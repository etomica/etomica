package etomica.action;

import java.util.Iterator;

import etomica.DataManager;
import etomica.Integrator;
import etomica.Phase;
import etomica.Simulation;

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
     * this method is broken
     */
    public void actionPerformed() {
        
        for(Iterator iter=simulation.getIntegratorList().iterator(); iter.hasNext(); ) {
            Integrator integrator = (Integrator)iter.next();
            integrator.initialize();
        }
        
        for(Iterator iter=simulation.getDataManagerList().iterator(); iter.hasNext(); ) {
            DataManager dataManager = (DataManager)iter.next();
            dataManager.resetAccumulators();
        }
    }
    
}
        