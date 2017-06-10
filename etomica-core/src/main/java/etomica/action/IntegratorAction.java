/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.integrator.Integrator;

/**
  * Elementary action performed on an integrator.
  */
public interface IntegratorAction extends IAction {

    public void setIntegrator(Integrator integrator);
    
    public Integrator getIntegrator();

}
