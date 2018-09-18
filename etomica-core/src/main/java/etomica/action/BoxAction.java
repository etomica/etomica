/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.box.Box;

/**
 * Elementary action performed on a Box. If action is invoked without previously specifying Box, no action is performed.
 */
public interface BoxAction extends IAction {

    /**
     * @return the Box specified by the last call to setBox, or null if no Box was previously specified
     */
    Box getBox();

    /**
     * Specifies the Box that will be subject to the action, upon subsequent call to actionPerformed()
     * @param box the Box targeted by the action
     */
    void setBox(Box box);

}
