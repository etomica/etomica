/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.atom.iterator.IteratorDirective.Direction;

/**
 * Adapater class that wraps three atomPair iterators, one suitable for
 * iterating over all molecule pairs in a box (AA), another suitable for
 * iterating over all molecule pairs formed with a target molecule (1A), and the
 * third suitable for iterating over a single molecule pair (11). Appropriate
 * iterator is selected based on argument given to setTarget method. If
 * setTarget method is never called, or has count == 0, the all-pair iterator is
 * used by default. <p>
 * This class may be set up to do inter- or intra-species iteration, depending
 * on choice of inner iterators given at construction.  Wrapped iterators are final
 * and cannot be changed after construction.
 * 
 */
public class MpiMolecule implements MoleculesetIteratorPDT, java.io.Serializable {

    /**
     * Constructs iterator by wrapping three others.
     * @param api1A
     *            iterator for all pairs formed with a target molecule
     * @param apiAA
     *            iterator for all pairs in the box
     */
    public MpiMolecule(MoleculesetIteratorPDT api1A, MoleculesetIteratorBoxDependent apiAA) {
        this.api1A = api1A;
        this.apiAA = apiAA;
        setTarget(null);
    }

    /**
     * Sets target atoms and determines the pair iterator according to the value
     * of targetAtom. Thus, for targetAtom equal to
     * <ul>
     * <li>null, AA iterator is indicated
     * <li>not-null, 1A iterator is indicated
     * </ul>
     * Target is passed on to setTarget method of corresponding iterator, which
     * ultimately determines the behavior of this iterator.
     */
    public void setTarget(IMolecule targetAtom) {
        if (targetAtom == null) {
            iterator = apiAA;
        }
        else {
            iterator = api1A;
            api1A.setTarget(targetAtom);
        }
    }

    /**
     * Specifies the box from which iterates are taken.
     */
    public void setBox(Box box) {
        this.box = box;
        if (box == null) {
            throw new NullPointerException("Null Box");
        }
    }

    /**
     * Specifies the direction from which partners of a target atoms
     * are taken to form iterates.  Has no effect if performing AA or 11
     * iteration (that is, if set target has count not equal to 1).
     */
    public void setDirection(Direction direction) {
        api1A.setDirection(direction);
    }

    /**
     * Readies iterator for iteration.
     */
    public void reset() {
        iterator.setBox(box);
        iterator.reset();
    }

    /**
     * Sets iterator to state in which hasNext is false.
     */
    public void unset() {
        apiAA.unset();
        api1A.unset();
    }

    /**
     * Same as nextPair.
     */
    public IMoleculeList next() {
        return iterator.next();
    }

    /**
     * Returns the number of iterates that would be given
     * when iterating after reset.  Does not require reset;
     * clobbers iteration state.
     */
    public int size() {
        iterator.setBox(box);
        return iterator.size();
    }

    /**
     * Returns 2, indicating that this is a pair iterator.
     */
    public final int nBody() {
        return 2;
    }

    /**
     * Returns the 1A iterator set at construction.
     */
    public MoleculesetIteratorPDT getApi1A() {
        return api1A;
    }

    /**
     * Returns the AA iterator set at construction.
     */
    public MoleculesetIteratorBoxDependent getApiAA() {
        return apiAA;
    }
    
    private static final long serialVersionUID = 1L;
    private MoleculesetIteratorBoxDependent iterator;
    private final MoleculesetIteratorPDT api1A;
    private final MoleculesetIteratorBoxDependent apiAA;
    private Box box;

}
