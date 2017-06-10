/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.AtomPair;

/**
 * Pair iterator synthesized from two atom iterators, such that the inner-loop
 * iteration depends on the outer-loop atom. Pairs are formed from the atoms
 * yielded by the two atom iterators. The inner-loop iterator must implement
 * AtomIteratorAtomDependent, and its set(Atom) method will be invoked with the
 * current outer-loop atom at the start of each inner-loop iteration. All pairs
 * returned by iterator are the same Atom[] instance, and differ only in the
 * Atom instances held by it.
 * <p>
 * Iterator can be condition to put the atoms in either order in the AtomPair that
 * it returns.  Thus the inner-loop Atom may be atom0 of the returned AtomPair,
 * or it may be atom1.  This behavior is set at construction, and cannot be changed
 * afterwards.  Default behavior has outer loop atoms as atom0, and inner loop atoms
 * as atom1.
 */
public class ApiInnerVariable implements AtomLeafsetIterator, java.io.Serializable {

    /**
     * Construct a pair iterator using the given atom iterators. Requires call
     * to reset() before beginning iteration.
     */
    public ApiInnerVariable(AtomIterator aiOuter,
            AtomIteratorAtomDependent aiInner) {
        this(aiOuter,aiInner,false);
    }
    
    /**
     * Construct a pair iterator using the given atom iterators, indicating
     * whether the atoms ordering in the AtomPair should be swapped from the
     * default behavior.
     * 
     * @param aiOuter
     *            outer-loop iterator
     * @param aiInner
     *            inner-loop iterator
     * @param doSwap
     *            if false (default), outer-loop atoms are given in atom0, and
     *            inner loop in atom1; if true outer-loop atoms are given in
     *            atom1, and inner loop in atom0
     */
    public ApiInnerVariable(AtomIterator aiOuter,
            AtomIteratorAtomDependent aiInner, boolean doSwap) {
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
    public AtomIterator getOuterIterator() {
        return aiOuter;
    }

    /**
     * Accessor method for the inner-loop atom iterator.
     * 
     * @return the current inner-loop iterator
     */
    public AtomIterator getInnerIterator() {
        return aiInner;
    }

    /**
     * Sets the iterator such that hasNext is false.
     */
    public void unset() {
        aiOuter.unset();
        aiInner.unset();
    }

    /**
     * Returns the number of pairs given by this iterator. Independent of state
     * of hasNext. Clobbers the iteration state (i.e., status of hasNext/next)
     * but does not recondition iterator (i.e., does not change set of iterates
     * that would be given on iteration after reset). Must perform reset if
     * attempting iteration after using size() method.
     */
    public int size() {
        int sum = 0;
        aiOuter.reset();
        for (IAtom a = aiOuter.nextAtom(); a != null; a = aiOuter.nextAtom()) {
            aiInner.setAtom(a);
            sum += aiInner.size();
        }
        return sum;
    }

    /**
     * Resets the iterator, so that it is ready to go through all of its pairs.
     */
    public void reset() {
        aiOuter.reset();
        IAtom nextOuter = aiOuter.nextAtom();
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
        aiInner.setAtom(nextOuter);
        aiInner.reset();
    }

    /**
     * Returns the next pair of atoms. The same Atom[] instance is returned
     * every time, but the Atoms it holds are (of course) different for each
     * iterate.
     */
    public IAtomList next() {
        IAtom nextInner = aiInner.nextAtom();
        while (nextInner == null) {
            IAtom nextOuter = aiOuter.nextAtom();
            if (nextOuter == null) {
                return null;
            }
            if (doSwap) {
                pair.atom1 = nextOuter;
            }
            else {
                pair.atom0 = nextOuter;
            }
            aiInner.setAtom(nextOuter);
            aiInner.reset();
            nextInner = aiInner.nextAtom();
            if (nextInner == null) {
                return null;
            }
        }
        
        if (doSwap) {
            pair.atom0 = nextInner;
        }
        else {
            pair.atom1 = nextInner;
        }
        
        if (pair.atom0 == null || pair.atom1 == null) {
            throw new RuntimeException("oops "+pair);
        }
        
        return pair;
    }

    public final int nBody() {
        return 2;
    }

    private static final long serialVersionUID = 2L;
    protected final AtomPair pair = new AtomPair();

    /**
     * The iterators used to generate the sets of atoms.
     */
    protected final AtomIteratorAtomDependent aiInner;
    protected final AtomIterator aiOuter;
    protected final boolean doSwap;
}

