/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.IController;
import etomica.integrator.Integrator;
import etomica.simulation.Simulation;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.exception.ConfigurationOverlapException;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeCubicSimple;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.space.Space;

/**
 * Action that invokes reset method of all registered simulation elements,
 * effectively initializing the entire simulation.
 */
public final class SimulationRestart extends SimulationActionAdapter {
    
    public SimulationRestart(Simulation sim, Space _space, IController _controller) {
        setSimulation(sim, _space, _controller);
    }

    protected void setSimulation(Simulation sim, Space _space, IController _controller) {
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

        IAction myAction = null;
        IAction[] currentActions = controller.getCurrentActions();
        if (currentActions.length == 1) {
            myAction = currentActions[0];
        }
        else if (currentActions.length == 0) {
            // we've reset the controller, which turns all the "current" actions to "pending"
            IAction[] pendingActions = controller.getPendingActions();
            if (pendingActions.length == 1) {
                myAction = pendingActions[0];
            }
        }
        if (myAction instanceof ActivityIntegrate) {
            Integrator integrator = ((ActivityIntegrate)myAction).getIntegrator();
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
