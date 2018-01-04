/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list;

import java.util.Arrays;


/**
 * Class used to maintain list of whether each 1-body potential that applies to an AtomType is currently in effect or
 * not for an Atom.
 */
public final class AtomPotentialList {

    private static final boolean[] EMPTY_BOOLEAN_ARR = new boolean[0];
    private boolean[] isInteractingList;

    /**
     * Constructs sequencer for the given atom.
     */
    public AtomPotentialList() {
        // use static empty array to avoid initially allocating a separate empty array for each atom
        isInteractingList = EMPTY_BOOLEAN_ARR;
    }

    /**
     * Sets whether the atom is interacting or not with the potential corresponding to the given index.
     */
    public void setIsInteracting(boolean flag, int index) {
        if (!(index < isInteractingList.length)) {
            isInteractingList = Arrays.copyOf(isInteractingList, index + 1);
        }
        isInteractingList[index] = flag;
    }

    /**
     * Returns booleans indicating whether the Atom interacts with each potential or not.
     */
    public boolean[] getInteractingList() {
        return isInteractingList;
    }

    protected void setCapacity(int index) {
        isInteractingList = Arrays.copyOf(isInteractingList, index);
    }

    /**
     * Clears neighbor lists, resetting all flags to false (not interacting).
     */
    public void clearNbrs() {
        final int length = isInteractingList.length;
        for (int i = 0; i < length; i++) {
            isInteractingList[i] = false;
        }
    }
}
