/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.api.IBox;

 /**
  * Elementary action performed on a box.
  */
public interface BoxAction extends IAction {

    public void setBox(IBox box);
    
    public IBox getBox();

}