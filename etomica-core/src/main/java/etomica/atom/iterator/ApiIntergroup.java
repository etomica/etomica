/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.AtomPair;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeSetSinglet;
import etomica.util.Debug;

/**
 * Iterator that returns pairs formed using two different basis atoms, so that
 * the iterates are taken from two different groups.
 */
public class ApiIntergroup implements AtomsetIteratorBasisDependent {

    /**
     * Default iterator is an ApiInnerFixed formed from two AtomIteratorBasis
     * instances.
     */
    public ApiIntergroup() {
        this(new AtomIteratorBasis(), new AtomIteratorBasis());
    }

    /**
     * Constructs a pair iterator that returns iterates from the given
     * pairIterator, which is expected to contain two basis-dependent 
     * iterators.
     */
    public ApiIntergroup(AtomIteratorBasisDependent outer, AtomIteratorBasisDependent inner) {
        super();
        aiOuter = outer;
        aiInner = inner;
        unset();
        atomSetSinglet = new MoleculeSetSinglet();
    }

    /**
     * Accessor method for the outer-loop atom iterator.
     * 
     * @return the current outer-loop iterator
     */
    public AtomIteratorBasisDependent getOuterIterator() {
        return aiOuter;
    }

    /**
     * Accessor method for the inner-loop atom iterator.
     * 
     * @return the current inner-loop iterator
     */
    public AtomIteratorBasisDependent getInnerIterator() {
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

    public final int nBody() {
        return 2;
    }

    /**
     * Returns the next pair of atoms. The same AtomPair instance is returned
     * every time, but the Atoms it holds are (of course) different for each
     * iterate.
     */
    public IAtomList next() {
        //Advance the inner loop, if it is not at its end.
        IAtom nextInner = aiInner.nextAtom();
        if (nextInner != null) {
            pair.atom1 = nextInner;
        }
        //Advance the outer loop, if the inner loop has reached its end.
        else {
            IAtom nextOuter = aiOuter.nextAtom();
            if (nextOuter == null) {
                return null;
            }
            aiInner.reset();
            nextInner = aiInner.nextAtom();
            if (nextInner == null) {
                return null;
            }

            pair.atom0 = nextOuter;
            pair.atom1 = nextInner;
        }
        if (Debug.ON && pair.atom0 == pair.atom1) {
            throw new RuntimeException("wow! "+pair);
        }
        return pair;
    }

    /**
     * Specifies a target atom, which should appear in all iterates. A
     * null value removes any restriction on the iterates.
     */
    public void setTarget(IAtom newTargetAtom) {
        targetAtom = newTargetAtom;
        needSetupIterators = true;
    }

    public boolean haveTarget(IAtom target) {
        if (target == null) {
            return true;
        }
        return aiOuter.haveTarget(target) || aiInner.haveTarget(target);
    }

    protected void setupIterators() {
        if (targetAtom == null) {
            aiOuter.setTarget(targetAtom);
            aiInner.setTarget(targetAtom);
        }
        else {
            if (aiInner.haveTarget(targetAtom)) {
                aiOuter.setTarget(null);
                aiInner.setTarget(targetAtom);
            } else {
                aiOuter.setTarget(targetAtom);
                aiInner.setTarget(null);
            }
        }
        needSetupIterators = false;
    }

    /**
     * Specifies the basis, which identifies the atoms subject to iteration. The
     * given atomSet must be of length 2. The first atom in the set specifies
     * the basis for the outer-loop iteration, and second atom specifies the
     * basis for the inner-loop iteration. In each case, if the basis atom is
     * not a leaf atom, its children will be the subject of iteration. If the
     * basis atom is a leaf, it will itself be the iterate. If given atomset is
     * null iterator will give no iterates until a proper basis is specified
     * via another call to this method.
     */
    public void setBasis(IMoleculeList basisAtoms) {
        if (basisAtoms == null) {
            aiOuter.setBasis(null);
        } else {
            atomSetSinglet.atom = basisAtoms.get(0);
            aiOuter.setBasis(atomSetSinglet);
            atomSetSinglet.atom = basisAtoms.get(1);
            aiInner.setBasis(atomSetSinglet);
        }
        needSetupIterators = true;
    }

    public void reset() {
        if (needSetupIterators) {
            setupIterators();
        }

        aiOuter.reset();
        aiInner.reset();
        IAtom nextOuter = aiOuter.nextAtom();
        if (nextOuter == null) {
            aiInner.unset();
            return;
        }

        pair.atom0 = nextOuter;
    }

    /**
     * Returns 2, indicating that the setBasis method expects an array of two
     * atoms.
     */
    public int basisSize() {
        return 2;
    }

    protected final AtomPair pair = new AtomPair();
    protected final AtomIteratorBasisDependent aiInner, aiOuter;
    protected IAtom targetAtom;
    protected boolean needSetupIterators = true;
    protected final MoleculeSetSinglet atomSetSinglet;

}
