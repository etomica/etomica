package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.utility.Arrays;

/**
 * Iterates over all the atoms given across an array of iterators. Atoms given
 * by multiple iterators, or by iterators added multiple times to this iterator,
 * will be iterated multiple times.
 * 
 * @author David Kofke
 */

public class AtomIteratorCompound implements AtomIterator, java.io.Serializable {

    /**
     * Construct iterator to loop over all iterates obtained from each iterator
     * in the given array. Must reset before using. Given array is cloned for
     * internal use, but references to iterators it contains are maintained
     * (iterators are not copied).
     */
    public AtomIteratorCompound(AtomIterator[] iterators) {
        iteratorSet = (AtomIterator[]) iterators.clone();
        unset();
    }

    /**
     * Constructs iterator that will give no iterates until an iterator is
     * added.
     */
    public AtomIteratorCompound() {
        iteratorSet = new AtomIterator[0];
        unset();
    }

    /**
     * Adds the given iterator to the set of iterators collected by this
     * compound iterator. No check is made the iterator was not already added;
     * if an iterator is added more than once, it will be iterated for each
     * addition, as if it were a different instance each time. This iterator is
     * left in an unset state.
     */
    public void addIterator(AtomIterator iterator) {
        if (iterator == null) {
            return;
        }
        iteratorSet = (AtomIterator[]) Arrays.addObject(iteratorSet, iterator);
        unset();
    }

    /**
     * Removes the given iterator from the set of iterators collected by this
     * compound iterator. If iterator was not previously added (or is null), no
     * action is taken. If iterator was added multiple times, only one reference
     * to it is removed. This iterator is left in an unset state.
     */
    public void removeIterator(AtomIterator iterator) {
        iteratorSet = (AtomIterator[]) Arrays.removeObject(iteratorSet,
                iterator);
        unset();
    }

    /**
     * Sets the iterators that will be iterated, discarding any iterators
     * previously added. If argument is null, no iterates will be given.
     * Reference to iterators is maintained, so changes to them will be
     * reflected by this iterator; reference to array is not kept (clone of
     * array is used internally by this iterator).
     */
    public void setIterators(AtomIterator[] iterators) {
        if (iterators == null) {
            iteratorSet = new AtomIterator[0];
        } else {
            iteratorSet = (AtomIterator[]) iterators.clone();
        }
        unset();
    }

    /**
     * Indicates whether iterator has another iterate.
     */
    public boolean hasNext() {
        return hasNext;
    }

    /**
     * Puts iterator in state in which hasNext is false.
     */
    public void unset() {
        hasNext = false;
    }

    /**
     * Returns 1, indicating that this is an atom iterator.
     */
    public int nBody() {
        return 1;
    }

    /**
     * Returns the number of iterates that would be given on reset of iterator.
     * Equal to the sum of the sizes of all iterators added to this compound
     * iterator.
     */
    public int size() {
        if (iteratorSet == null) {
            return 0;
        }
        int count = 0;
        for (int i = 0; i < iteratorSet.length; i++) {
            count += iteratorSet[i].size();
        }
        return count;
    }

    /**
     * Returns true if the given atom is among the iterates given by this
     * iterator. Does not affect iteration state. Returns false if argument is
     * null or if atoms.count != 1.
     */
    public boolean contains(AtomSet atoms) {
        if (iteratorSet == null) {
            return false;
        }
        for (int i = 0; i < iteratorSet.length; i++) {
            if (iteratorSet[i].contains(atoms)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Readies iterator to begin iteration.
     */
    public void reset() {
        if (!hasIterator()) {
            hasNext = false;
            return;
        }
        index = 0;
        iteratorSet[index].reset();
        while (!iteratorSet[index].hasNext() && index + 1 < iteratorSet.length) {
            index++;
            iteratorSet[index].reset();
        }
        hasNext = iteratorSet[index].hasNext();
    }

    /**
     * Returns next iterate without advancing the iterator.
     */
    public AtomSet peek() {
        return hasNext ? iteratorSet[index].peek() : null;
    }

    /**
     * Same as nextAtom.
     */
    public AtomSet next() {
        return nextAtom();
    }

    /**
     * Returns next iterate and advances iterator.
     */
    public Atom nextAtom() {
        if (!hasNext) {
            return null;
        }
        Atom atom = iteratorSet[index].nextAtom();
        while (!iteratorSet[index].hasNext()) {
            if (++index < iteratorSet.length) {
                iteratorSet[index].reset();
            } else {
                hasNext = false;
                break;
            }
        }
        return atom;
    }

    /**
     * Performs action on all atoms in added iterators. Performs action multiple
     * times to atoms present in multiple iterators, or given by iterators added
     * multiple times to compound iterator.
     */
    public void allAtoms(AtomsetAction action) {
        for (int i = 0; i < iteratorSet.length; i++) {
            iteratorSet[i].allAtoms(action);
        }
    }

    private boolean hasIterator() {
        return (iteratorSet.length > 0);
    }

    private AtomIterator[] iteratorSet = new AtomIterator[0];
    private boolean hasNext;
    private int index;

}//end of AtomIteratorCompound
