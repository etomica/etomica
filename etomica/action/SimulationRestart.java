package etomica.action;

import java.util.Iterator;

import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.config.ConfigurationSequential;
import etomica.data.DataAccumulator;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.Integrator;
import etomica.lattice.LatticeCubicFcc;
import etomica.phase.Phase;
import etomica.simulation.Simulation;

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
        if (sim.space().D() == 3) {
            setConfiguration(new ConfigurationLattice(new LatticeCubicFcc()));
        }
        else {
            setConfiguration(new ConfigurationSequential(sim.space()));
        }
        ignoreOverlap = sim.getDefaults().ignoreOverlap;
        accumulatorAction = new SimulationDataAction(sim, new ResetAccumulators());
    }
        
    /**
     * Resets phases, integrators, and accumulators.
     */
    public void actionPerformed() {

        Phase[] phases = simulation.getPhases();
        for(int i=0; i<phases.length; i++) {
            if (configuration != null) {
                configuration.initializeCoordinates(phases[i]);
            }
        }
        
        for(Iterator iter=simulation.getIntegratorList().iterator(); iter.hasNext(); ) {
            Integrator integrator = (Integrator)iter.next();
            if(integrator.isInitialized()) {
                try {
                    integrator.initialize();
                }
                catch (ConfigurationOverlapException e) {
                    if (!ignoreOverlap) {
                        throw new RuntimeException(e);
                    }
                }
            }
        }
    
        accumulatorAction.actionPerformed();
    }
    
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

    private Configuration configuration;
    private boolean ignoreOverlap;
    private final SimulationDataAction accumulatorAction;
}