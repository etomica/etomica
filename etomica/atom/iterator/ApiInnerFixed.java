package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.atom.Atom;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;

/**
 * Pair iterator synthesized from two atom iterators, such that the inner-loop
 * iteration is independent of the outer-loop atom. Pairs are formed from the
 * atoms yielded by the two atom iterators. It is expected that the inner-loop
 * iterator will yield the same set of atoms with each pass of the outer loop.
 * All pairs returned by iterator are the same AtomPair instance, and differ only
 * in the Atom instances held by it.
 */

/*
 * History of changes 08/25/04 (DAK et al) new
 */

public final class ApiInnerFixed implements AtomPairIterator, ApiComposite, java.io.Serializable {

    /**
     * Construct a pair iterator using the given atom iterators. Requires call
     * to reset() before beginning iteration.
     */
    public ApiInnerFixed(AtomIterator aiOuter, AtomIterator aiInner) {
        this.aiOuter = aiOuter;
        this.aiInner = aiInner;
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
        hasNext = false;
    }

    /**
     * Indicates whether the given atom pair will be returned by the iterator
     * during its iteration. The order of the atoms in the pair is significant
     * (this means that a value of true is returned only if one of the pairs
     * returned by the iterator will have the same two atoms in the same
     * atom1/atom2 position as the given pair). Not dependent on state of
     * hasNext. Returns false if pair is null.
     */
    public boolean contains(AtomSet pair) {
        if (pair == null || pair.count() != 2) {
            return false;
        }
        return aiOuter.contains(pair.getAtom(0))
                && aiInner.contains(pair.getAtom(1));
    }

    /**
     * Returns the number of pairs given by this iterator. Not dependent on
     * state of hasNext.
     */
    public int size() {
        return aiOuter.size() * aiInner.size();
    }

    /**
     * Indicates whether the iterator has completed its iteration.
     */
    public boolean hasNext() {
        return hasNext;
    }

    /**
     * Resets the iterator, so that it is ready to go through all of its pairs.
     * A previously returned pair may be altered by this method.
     */
    public void reset() {
        aiOuter.reset();
        aiInner.reset();
        hasNext = aiOuter.hasNext() && aiInner.hasNext();
        if (hasNext) {
            pair.atom0 = aiOuter.nextAtom();
        }
    }

    /**
     * Returns the next pair without advancing the iterator. If the iterator has
     * reached the end of its iteration, returns null. A previously-returned
     * pair will be altered by this method.
     */
    public AtomSet peek() {
        if (!hasNext) {
            return null;
        }

        if (aiInner.hasNext()) {
            pair.atom1 = (Atom) aiInner.peek();
        } else {
            // Althouth we advance aiOuter, we
            // are not advancing the pair iterator.
            // Outcome of next() is not changed
            aiInner.reset();
            pair.atom0 = aiOuter.nextAtom();
            pair.atom1 = (Atom) aiInner.peek();
        }
        return pair;
    }

    public final AtomSet next() {
        return nextPair();
    }

    /**
     * Returns the next pair of atoms. The same AtomPair instance is returned
     * every time, but the Atoms it holds are (of course) different for each
     * iterate.
     */
    public AtomPair nextPair() {
        if (!hasNext) {
            return null;
        }
        //Advance the inner loop, if it is not at its end.
        if (aiInner.hasNext()) {
            pair.atom1 = aiInner.nextAtom();
        }
        //Advance the outer loop, if the inner loop has reached its end.
        else {
            aiInner.reset();
            pair.atom0 = aiOuter.nextAtom();
            pair.atom1 = aiInner.nextAtom();
        }
        hasNext = aiInner.hasNext() || aiOuter.hasNext();
        return pair;
    }

    /**
     * Performs the given action on all pairs returned by this iterator.
     */
    public void allAtoms(AtomsetAction action) {
        aiOuter.reset();
        while (aiOuter.hasNext()) {
            pair.atom0 = aiOuter.nextAtom();
            aiInner.reset();
            while (aiInner.hasNext()) {
                pair.atom1 = aiInner.nextAtom();
                action.actionPerformed(pair);
            }
        }
    }

    public final int nBody() {
        return 2;
    }

    private final AtomPair pair = new AtomPair();
    private boolean hasNext;
    private final AtomIterator aiInner, aiOuter;

} //end of class ApiInnerFixed

