package etomica.action;

import etomica.*;
import java.util.Iterator;

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
        
        simulation.getController().halt();
        try {
            simulation.getController().wait();
        }
        catch (InterruptedException e) {}
        
        for(Iterator iter=simulation.getPhaseList().iterator(); iter.hasNext(); ) {
            Phase phase = (Phase)iter.next();
            phase.getConfiguration().initializeCoordinates((phase.speciesMaster().node).childAtomArray());
        }
        
        for(Iterator iter=simulation.getIntegratorList().iterator(); iter.hasNext(); ) {
            Integrator integrator = (Integrator)iter.next();
            integrator.reset();
        }
        
        for(Iterator iter=simulation.getDataManagerList().iterator(); iter.hasNext(); ) {
            DataManager dataManager = (DataManager)iter.next();
            dataManager.resetAccumulators();
        }
        
/*        for(Iterator iter=simulation.getDisplayList().iterator(); iter.hasNext(); ) {
            Display display = (Display)iter.next();
            display.doUpdate();
            display.graphic().repaint();
        }*/
    }
    
}
        