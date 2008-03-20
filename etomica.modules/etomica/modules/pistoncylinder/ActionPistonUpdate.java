package etomica.modules.pistoncylinder;

import etomica.api.IAction;

/**
 * Action that informs the integrator that some property of the piston
 * has changed (pressure, mass, velocity)
 */

public class ActionPistonUpdate implements IAction {
    public ActionPistonUpdate(IntegratorHardPiston integrator) {
        pistonIntegrator = integrator;
    }
    
    public void actionPerformed() {
        pistonIntegrator.pistonUpdateRequested();
    }
    
    public String getLabel() {return "Piston updater";}
    
    private final IntegratorHardPiston pistonIntegrator;
}