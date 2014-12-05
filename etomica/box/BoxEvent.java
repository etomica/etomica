/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.api.IBox;
import etomica.api.IBoxEvent;

/**
 * Event that conveys some happening with respect to a box or the things it contains.
 *
 * @see BoxListenerAdapter
 */
public class BoxEvent implements java.io.Serializable, IBoxEvent {
    
    public BoxEvent(IBox box) {
        this.box = box;
    }
    
    public IBox getBox() {
        return box;
    }
    
    protected IBox box;

    private static final long serialVersionUID = 1L;
}
    