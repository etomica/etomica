package etomica.action;

import java.util.Iterator;

import etomica.Simulation;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.config.ConfigurationSequential;
import etomica.data.DataAccumulator;
import etomica.integrator.Integrator;
import etomica.lattice.LatticeCubicFcc;
import etomica.phase.Phase;

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
        if (sim.space().D() == 2) {
            setConfiguration(new ConfigurationSequential(sim.space()));
        }
        else {
            setConfiguration(new ConfigurationLattice(new LatticeCubicFcc()));
        }
    }
        
    /**
     * Resets phases, integrators, and accumulators.
     */
    public void actionPerformed() {

        for(Iterator iter=simulation.getPhaseList().iterator(); iter.hasNext(); ) {
            if (configuration != null) {
                configuration.initializeCoordinates((Phase)iter.next());
            }
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
    
    private Configuration configuration;
    /**
     * @return Returns the configuration.
     */
    public Configuration getConfiguration() {
        return configuration;
    }
    /**
     * @param configuration The configuration to set.
     */
    public void setConfiguration(Configuration configuration) {
        this.configuration = configuration;
    }
}
        