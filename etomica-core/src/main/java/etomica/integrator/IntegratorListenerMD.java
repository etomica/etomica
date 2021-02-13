/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

public interface IntegratorListenerMD extends IntegratorListener {

    /**
     * Invoked before the integrator computes the forces on all atoms
     * in the system.
     *
     * @param e
     */
    default void integratorForcePrecomputed(IntegratorEvent e) {
    }

    /**
     * Invoked after the integrator has computed the forces on all atoms
     * in the system.
     *
     * @param e
     */
    default void integratorForceComputed(IntegratorEvent e) {
    }

}
