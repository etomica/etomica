/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;

import etomica.box.Box;

/**
 * A box event fired by the box's event manager.  The box that changed is
 * available via getBox.  Other details can be inferred from the other
 * interfaces implemented by the event object and/or obtained from calling
 * methods on those interfaces.
 */
public interface IBoxEvent {
    
    /**
     * @return the box that changed and triggered the event
     */
    public Box getBox();

}
