package etomica.action;

import java.util.Iterator;

import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.Integrator;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeCubicSimple;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.phase.Phase;
import etomica.simulation.Simulation;

/**
 * Action that invokes reset method of all registered simulation elements,
 * effectively initializing the entire simulation.
 */
public final class SimulationRestart extends SimulationActionAdapter {
    
    public SimulationRestart(Simulation sim) {
        setSimulation(sim);
    }
    
    public void setSimulation(Simulation sim) {
        super.setSimulation(sim);
        if (sim.getSpace().D() == 3) {
            setConfiguration(new ConfigurationLattice(new LatticeCubicFcc()));
        }
        else if (sim.getSpace().D() == 2) {
            setConfiguration(new ConfigurationLattice(new LatticeOrthorhombicHexagonal()));
        }
        else {
            setConfiguration(new ConfigurationLattice(new LatticeCubicSimple(1, 1.0)));
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

    private static final long serialVersionUID = 1L;
    private Configuration configuration;
    private boolean ignoreOverlap;
    private SimulationDataAction accumulatorAction;
}