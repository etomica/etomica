/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.api.IBox;


/**
 * Event that conveys that the maximum global index in a Box has changed.
 */
public class BoxGlobalAtomIndexEvent extends BoxEvent {

    public BoxGlobalAtomIndexEvent(IBox box, int maxIndex) {
        super(box);
        this.maxIndex = maxIndex;
    }
    
    public int getMaxIndex() {
        return maxIndex;
    }
    
    private final int maxIndex;
    private static final long serialVersionUID = 1L;
}
