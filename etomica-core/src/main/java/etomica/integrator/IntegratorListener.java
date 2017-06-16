/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

public interface IntegratorListener {

    /**
     * Invoked when integration begins.
     *
     * @param e
     */
    void integratorInitialized(IntegratorEvent e);

    /**
     * Invoked at the beginning of each integrator step.
     *
     * @param e
     */
    void integratorStepStarted(IntegratorEvent e);

    /**
     * Invoked at the end of each integrator step.
     *
     * @param e
     */
    void integratorStepFinished(IntegratorEvent e);

}
