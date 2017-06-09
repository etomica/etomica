/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.*;
import etomica.box.Box;
import etomica.potential.IteratorDirective.Direction;

/**
 * Iterator for all the molecules of a set of species in a box.  Each iterate
 * is all the molecules in a box, with each Atom as the first atom in the 
 * set. This class is used by PotentialMaster to iterate over molecules for 
 * N-body potentials.
 * 
 * This class is designed to work and conform to the API... not to be efficient 
 * or pleasant to look at!  Use neighbor lists. 
 */
public class AtomIteratorAllLeafType implements AtomsetIteratorBoxDependent,
                AtomsetIteratorDirectable, AtomsetIteratorTargetable, java.io.Serializable {

    private static final long serialVersionUID = 1L;
    private final AtomType[] atomType;
    private final AtomListWrapper next;
    private Box box;
    private int nextCursor;
    
    /**
     * @param atomType species for which molecules are returned as iterates. Only
     * species[0] is relevant, and must not be null.
     */
    public AtomIteratorAllLeafType(AtomType[] atomType) {
        this.atomType = atomType;
        next = new AtomListWrapper();
    }

    /**
     * Sets the box containing the molecules for iteration. A null
     * box conditions iterator to give no iterates.
     */
    public void setBox(Box newBox) {
        box = newBox;
        if (box == null) {
            throw new NullPointerException("Null box");
        }
    }

    /**
     * Sets the target of iteration... has no actual effect since all iterates
     * contain all Atoms.
     */
    public void setTarget(IAtom newTargetAtom) {
    }

    /**
     * Has no effect, but is included as part of the AtomsetIteratorPDT interface.
     * Besides, you didn't really want to iterate down, did you?
     */
    public void setDirection(Direction newDirection) {
    }

    public void reset() {
        // add all Atoms to ArrayList we will return
        AtomArrayList arrayList = next.getArrayList();
        arrayList.clear();
        IAtomList leafList = box.getLeafList();
        for (int i=0; i<leafList.getAtomCount(); i++) {
        	for (int j=0; j<atomType.length; j++) {
        		if(leafList.getAtom(i).getType()==atomType[j]){
        			arrayList.add(leafList.getAtom(i));
        		}
        	}
        }
        nextCursor = 0;
    }

    public void unset() {
        next.getArrayList().clear();
    }

    public IAtomList next() {
        if (nextCursor + 1 > next.getAtomCount()) {
            return null;
        }
        if (nextCursor < 0) {
            // already poked
            nextCursor = -nextCursor;
            return next;
        }
        AtomArrayList arrayList = next.getArrayList();
        IAtom oldFirst = arrayList.getAtom(0);
        arrayList.set(0,arrayList.getAtom(nextCursor));
        arrayList.set(nextCursor,oldFirst);
        nextCursor++;
        return next;
    }

    public int nBody() {
        return Integer.MAX_VALUE;
    }

    /**
     * Returns the number of iterates given by this iterator, if iterated after
     * a call to reset().
     */
    public int size() {
        return next.getAtomCount();
    }
}
