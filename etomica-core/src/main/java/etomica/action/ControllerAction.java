/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.action.activity.Controller;

/**
  * Elementary action performed on a controller.
  */
public interface ControllerAction extends IAction {

    public void setController(Controller c);
    
    public Controller getController();

}
