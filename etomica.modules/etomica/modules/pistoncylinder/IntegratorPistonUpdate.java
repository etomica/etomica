package etomica.modules.pistoncylinder;

import etomica.Action;

/**
 * Action that informs the integrator that some property of the piston
 * has changed (pressure, mass, velocity)
 */

class IntegratorPistonUpdate implements Action {
    public IntegratorPistonUpdate(IntegratorHardPiston integrator) {
        pistonIntegrator = integrator;
    }
    
    public void actionPerformed() {
        pistonIntegrator.pistonUpdateRequested();
    }
    
    public String getLabel() {return "Piston updater";}
    
    private final IntegratorHardPiston pistonIntegrator;
}