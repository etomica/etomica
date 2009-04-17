package etomica.action;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAction;
import etomica.api.IController;
import etomica.api.IIntegrator;
import etomica.api.ISimulation;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.exception.ConfigurationOverlapException;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeCubicSimple;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.space.ISpace;
import etomica.space.Space;

/**
 * Action that invokes reset method of all registered simulation elements,
 * effectively initializing the entire simulation.
 */
public final class SimulationRestart extends SimulationActionAdapter {
    
    public SimulationRestart(ISimulation sim, ISpace _space, IController _controller) {
        setSimulation(sim, _space, _controller);
    }

    protected void setSimulation(ISimulation sim, ISpace _space, IController _controller) {
        super.setSimulation(sim, _space);
        controller = _controller;
        if (space != null) {
            if (space.D() == 3) {
                setConfiguration(new ConfigurationLattice(new LatticeCubicFcc(space), space));
            }
            else if (space.D() == 2) {
                setConfiguration(new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space));
            }
            else {
            	Space sp = Space.getInstance(1);
                setConfiguration(new ConfigurationLattice(new LatticeCubicSimple(sp, 1.0), sp));
            }
        }
        ignoreOverlap = false;
        accumulatorAction = new SimulationDataAction(new ResetAccumulatorsAveraged());
    }

    public SimulationDataAction getDataResetAction() {
        return accumulatorAction;
    }
    
    public void setDataResetAction(SimulationDataAction newResetAction) {
        accumulatorAction = newResetAction;
    }

    public void setIgnoreOverlap(boolean doIgnoreOverlap) {
        ignoreOverlap = doIgnoreOverlap;
    }
    
    public boolean isIgnoreOverlap() {
        return ignoreOverlap;
    }
    
    /**
     * Resets boxs, integrators, and accumulators.
     */
    public void actionPerformed() {
        int boxCount = simulation.getBoxCount();
        for(int i=0; i<boxCount; i++) {
            if (configuration != null) {
                configuration.initializeCoordinates(simulation.getBox(i));
            }
        }

        IAction[] currentActions = controller.getCurrentActions();
        if (currentActions.length == 1) {
            IAction currentAction = currentActions[0];
            if (currentAction instanceof ActivityIntegrate) {
                IIntegrator integrator = ((ActivityIntegrate)currentAction).getIntegrator();
                if(integrator.getStepCount() > 0) {
                    integrator.resetStepCount();
                    try {
                        integrator.reset();
                    }
                    catch (ConfigurationOverlapException e) {
                        if (!ignoreOverlap) {
                            throw e;
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
    private IController controller;
}