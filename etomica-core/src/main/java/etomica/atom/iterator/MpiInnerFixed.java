/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.atom.MoleculePair;

/**
 * Pair iterator synthesized from two atom iterators, such that the inner-loop
 * iteration is independent of the outer-loop atom. Pairs are formed from the
 * atoms yielded by the two atom iterators. It is expected that the inner-loop
 * iterator will yield the same set of atoms with each pass of the outer loop.
 * All pairs returned by iterator are the same AtomPair instance, and differ only
 * in the Atom instances held by it.
 * <p>
 * Iterator can be condition to put the atoms in either order in the AtomPair that
 * it returns.  Thus the inner-loop Atom may be atom0 of the returned AtomPair,
 * or it may be atom1.  This behavior is set at construction, and cannot be changed
 * afterwards.  Default behavior has outer loop atoms as atom0, and inner loop atoms
 * as atom1.
 */
public final class MpiInnerFixed implements MoleculesetIterator, java.io.Serializable {

    /**
     * Construct a pair iterator using the given atom iterators. Requires call
     * to reset() before beginning iteration.
     */
    /**
     * Construct a pair iterator using the given atom iterators, indicating
     * whether the atoms ordering in the AtomPair should be swapped from the
     * default behavior.
     * 
     * @param aiOuter
     *            outer-loop iterator
     * @param aiInner
     *            inner-loop iterator
     */
    public MpiInnerFixed(MoleculeIterator aiOuter, MoleculeIterator aiInner,
            boolean doSwap) {
        this.aiOuter = aiOuter;
        this.aiInner = aiInner;
        this.doSwap = doSwap;
        unset();
    }

    /**
     * Accessor method for the outer-loop atom iterator.
     * 
     * @return the current outer-loop iterator
     */
    public MoleculeIterator getOuterIterator() {
        return aiOuter;
    }

    /**
     * Accessor method for the inner-loop atom iterator.
     * 
     * @return the current inner-loop iterator
     */
    public MoleculeIterator getInnerIterator() {
        return aiInner;
    }

    /**
     * Sets the iterator such that hasNext is false.
     */
    public void unset() {
        aiInner.unset();
        aiOuter.unset();
    }

    /**
     * Returns the number of pairs given by this iterator. Not dependent on
     * state of hasNext.
     */
    public int size() {
        return aiOuter.size() * aiInner.size();
    }

    /**
     * Resets the iterator, so that it is ready to go through all of its pairs.
     * A previously returned pair may be altered by this method.
     */
    public void reset() {
        aiOuter.reset();
        aiInner.reset();
        IMolecule nextOuter = aiOuter.nextMolecule();
        if (nextOuter == null) {
            aiInner.unset();
            return;
        }

        if (doSwap) {
            pair.atom1 = nextOuter;
        }
        else {
            pair.atom0 = nextOuter;
        }
    }

    /**
     * Returns the next pair of atoms. The same AtomPair instance is returned
     * every time, but the Atoms it holds are (of course) different for each
     * iterate.
     */
    public final IMoleculeList next() {
        //Advance the inner loop, if it is not at its end.
        IMolecule nextInner = aiInner.nextMolecule();
        if (nextInner != null) {
            if (doSwap) {
                pair.atom0 = nextInner;
            }
            else {
                pair.atom1 = nextInner;
            }
        }
        //Advance the outer loop, if the inner loop has reached its end.
        else {
            IMolecule nextOuter = aiOuter.nextMolecule();
            if (nextOuter == null) {
                return null;
            }
            aiInner.reset();
            nextInner = aiInner.nextMolecule();
            if (nextInner == null) {
                return null;
            }

            if (doSwap) {
                pair.atom1 = nextOuter;
                pair.atom0 = nextInner;
            }
            else {
                pair.atom0 = nextOuter;
                pair.atom1 = nextInner;
            }
        }
        if (pair.atom0 == pair.atom1) {
            throw new RuntimeException("wow! "+pair);
        }
        return pair;
    }

    public final int nBody() {
        return 2;
    }

    private static final long serialVersionUID = 2L;
    private final MoleculePair pair = new MoleculePair();
    private final MoleculeIterator aiInner, aiOuter;
    private final boolean doSwap;
}

