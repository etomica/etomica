/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

public class IntegratorEvent {

    private final Integrator integrator;

    public IntegratorEvent(Integrator integrator) {
        this.integrator = integrator;
    }

    public Integrator getIntegrator() {
        return integrator;
    }
}
