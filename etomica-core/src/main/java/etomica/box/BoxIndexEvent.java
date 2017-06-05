/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;


/**
 * Event that conveys that the maximum leaf index in a Box has changed
 * (or is about to change).
 */
public class BoxIndexEvent extends BoxEvent {

    protected final int index;

    public BoxIndexEvent(Box box, int index) {
        super(box);
        this.index = index;
    }

    public int getIndex() {
        return index;
    }
}
