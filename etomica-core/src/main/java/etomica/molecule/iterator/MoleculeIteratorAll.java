/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule.iterator;

import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;
import etomica.molecule.MoleculeListWrapper;
import etomica.potential.IteratorDirective.Direction;
import etomica.species.ISpecies;

/**
 * Iterator for all the molecules of a set of species in a box.  Each iterate
 * is all the molecules in a box, with each Atom as the first atom in the 
 * set. This class is used by PotentialMaster to iterate over molecules for 
 * N-body potentials.
 * 
 * This class is designed to work and conform to the API... not to be efficient 
 * or pleasant to look at!  Use neighbor lists. 
 */
public class MoleculeIteratorAll implements MoleculesetIteratorPDT, java.io.Serializable {

	public MoleculeIteratorAll(ISpecies[] species) {
		this(species,false);
	}
    /**
     * @param species species for which molecules are returned as iterates. Only
     * species[0] is relevant, and must not be null.
     */
    public MoleculeIteratorAll(ISpecies[] species,boolean oneIterate) {
        this.species = species;
        this.oneIterate=oneIterate;
        next = new MoleculeListWrapper();
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
    public void setTarget(IMolecule newTargetAtom) {
    }

    /** 
     * Has no effect, but is included as part of the AtomsetIteratorPDT interface.
     * Besides, you didn't really want to iterate down, did you?
     */
    public void setDirection(Direction newDirection) {
    }

    public void reset() {
        // add all Atoms to ArrayList we will return
        MoleculeArrayList arrayList = next.getArrayList();
        arrayList.clear();
        for (int i=0; i<species.length; i++) {
            arrayList.addAll(box.getMoleculeList(species[i]));
        }
        nextCursor = 0;
    }
    
    public void unset() {
        next.getArrayList().clear();
    }
    
    public IMoleculeList next() {
        if (nextCursor + 1 > next.size()||(oneIterate && nextCursor>0)) {
            return null;
        }
        if (nextCursor < 0) {
            // already poked
            nextCursor = -nextCursor;
            return next;
        }
        MoleculeArrayList arrayList = next.getArrayList();
        IMolecule oldFirst = arrayList.get(0);
        arrayList.set(0,arrayList.get(nextCursor));
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
        return next.size();
    }

    protected final ISpecies[] species;
    protected Box box;
    protected int nextCursor;
    protected final MoleculeListWrapper next;
    protected final boolean oneIterate;
}
