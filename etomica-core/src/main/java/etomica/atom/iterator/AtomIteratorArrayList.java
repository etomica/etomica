/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.*;
import etomica.potential.IteratorDirective;

/**
  * An atom iterator of the elements from an AtomArrayList (in proper
  * sequence).  Iterator will fail if element are added to or removed 
  * from list while iteration is proceeding.
  */
public class AtomIteratorArrayList extends AtomIteratorArrayListSimple implements AtomIteratorAtomDependent {

    /**
     * Constructs new iterator with an empty list.
     */
 	public AtomIteratorArrayList(IteratorDirective.Direction direction) {
 		this(direction, 0);
 	}
    
    /**
     * Constructs new iterator set to iterate given list (upon reset).
     */
 	public AtomIteratorArrayList(IteratorDirective.Direction direction, int numToSkip) {
        this(direction, numToSkip, new AtomToIndexChild(), new AtomToParentChildList());
    }
    
    public AtomIteratorArrayList(IteratorDirective.Direction direction, int numToSkip, 
            AtomToIndex atomToIndex, AtomToAtomLeafList atomToAtomSet) {
        if (direction == null)
            throw new IllegalArgumentException(
                    "Must specify direction to constructor of AtomLinkerIterator");
        upListNow = (direction == IteratorDirective.Direction.UP);

        this.atomToIndex = atomToIndex;
        this.atomToAtomSet = atomToAtomSet;
        this.numToSkip = numToSkip;
 	}
    
    /**
     * Returns the next iterate and advances the iterator.
     */
 	public IAtom nextAtom() {
        if (upListNow) {
            if (cursor > list.size()-2) {
                return null;
            }
            cursor++;
        }
        else {
            if (cursor < 1) {
                return null;
            }
            cursor--;
        }
 		return list.get(cursor);
 	}
 	
    /**
     * Returns the number of iterates that would be given by this iterator
     * if reset with the current list.
     */
    public int size() {
        int count = 0;
        reset();
        for (Object a = next(); a != null; a = next()) {
            count++;
        }
        return count;
    }

    /**
     * Puts iterator in state ready to begin iteration.
     */
 	public void reset() {
        list = atomToAtomSet.getAtomList(startAtom);
 		cursor = atomToIndex.getIndex(startAtom);
        if (upListNow) {
            cursor--;
            cursor += numToSkip;
        }
        else {
            cursor++;
            cursor -= numToSkip;
        }
 	}
    
    public void unset() {
        if (upListNow) {
            cursor = list.size();
        }
        else {
            cursor = -1;
        }
    }
 	
    /**
     * Sets the first atom for iteration. Iteration proceeds from this atom up
     * and/or down the list, depending on how iterator was configured at
     * construction.
     */
    public void setAtom(IAtom atom) {
        startAtom = atom;
    }

    private static final long serialVersionUID = 1L;
    protected final boolean upListNow;
    private final int numToSkip;
    private IAtom startAtom;
    private final AtomToIndex atomToIndex;
    private final AtomToAtomLeafList atomToAtomSet;
 }
