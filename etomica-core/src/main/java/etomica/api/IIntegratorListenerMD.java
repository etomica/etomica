/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;

import etomica.integrator.IntegratorEvent;

public interface IIntegratorListenerMD extends IIntegratorListener {
    
    /**
     * Invoked before the integrator computes the forces on all atoms
     * in the system.
     * @param e
     */
    public void integratorForcePrecomputed(IntegratorEvent e);
    
    /**
     * Invoked after the integrator has computed the forces on all atoms
     * in the system.
     * @param e
     */
    public void integratorForceComputed(IntegratorEvent e);
    
}
