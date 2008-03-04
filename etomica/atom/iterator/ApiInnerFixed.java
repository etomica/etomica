package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.atom.AtomPair;

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
public final class ApiInnerFixed implements ApiComposite, java.io.Serializable {

    /**
     * Construct a pair iterator using the given atom iterators. Requires call
     * to reset() before beginning iteration.
     */
    public ApiInnerFixed(AtomIterator aiOuter, AtomIterator aiInner) {
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
    public ApiInnerFixed(AtomIterator aiOuter, AtomIterator aiInner,
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
    }

    /**
     * Returns the next pair of atoms. The same AtomPair instance is returned
     * every time, but the Atoms it holds are (of course) different for each
     * iterate.
     */
    public final IAtomSet next() {
        //Advance the inner loop, if it is not at its end.
        IAtom nextInner = aiInner.nextAtom();
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
            IAtom nextOuter = aiOuter.nextAtom();
            if (nextOuter == null) {
                return null;
            }
            aiInner.reset();
            nextInner = aiInner.nextAtom();
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

    /**
     * Performs the given action on all pairs returned by this iterator.
     */
    public void allAtoms(AtomsetAction action) {
        aiOuter.reset();
        for (IAtom outerAtom = aiOuter.nextAtom(); outerAtom != null;
             outerAtom = aiOuter.nextAtom()) {
            if (doSwap) {
                pair.atom1 = outerAtom;
            }
            else {
                pair.atom0 = outerAtom;
            }
            aiInner.reset();
            for (IAtom innerAtom = aiInner.nextAtom(); innerAtom != null;
                 innerAtom = aiInner.nextAtom()) {
                if (doSwap) {
                    pair.atom0 = innerAtom;
                }
                else {
                    pair.atom1 = innerAtom;
                }
                action.actionPerformed(pair);
            }
        }
    }

    public final int nBody() {
        return 2;
    }

    private static final long serialVersionUID = 2L;
    private final AtomPair pair = new AtomPair();
    private final AtomIterator aiInner, aiOuter;
    private final boolean doSwap;
}

