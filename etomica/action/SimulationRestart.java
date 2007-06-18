package etomica.action;

import java.util.ArrayList;

import etomica.action.activity.ActivityIntegrate;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IIntegrator;
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
        integratorList = new ArrayList();
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
        ignoreOverlap = false;
        accumulatorAction = new SimulationDataAction(sim, new ResetAccumulators());
    }

    public void setIgnoreOverlap(boolean doIgnoreOverlap) {
        ignoreOverlap = doIgnoreOverlap;
    }
    
    public boolean isIgnoreOverlap() {
        return ignoreOverlap;
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
        
        Action[] currentActions = simulation.getController().getCurrentActions();
        if (currentActions.length == 1) {
            Action currentAction = currentActions[0];
            if (currentAction instanceof ActivityIntegrate) {
                IIntegrator integrator = ((ActivityIntegrate)currentAction).getIntegrator();
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
    protected Configuration configuration;
    protected boolean ignoreOverlap;
    protected SimulationDataAction accumulatorAction;
    protected ArrayList integratorList;
}