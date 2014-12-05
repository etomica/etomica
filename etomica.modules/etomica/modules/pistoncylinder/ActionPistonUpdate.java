/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.pistoncylinder;

import etomica.action.IAction;

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