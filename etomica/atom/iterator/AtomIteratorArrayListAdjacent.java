package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomList;
import etomica.atom.AtomSet;
import etomica.atom.AtomToArrayList;
import etomica.atom.AtomToIndex;
import etomica.atom.AtomToIndexChild;
import etomica.atom.AtomToParentChildList;
import etomica.atom.iterator.IteratorDirective.Direction;


/**
 * Returns one or both of the atoms adjacent to a specified atom in
 * its sequence list (i.e., the AtomList containing its seq linker).
 * If adjacent linker has no atom (is a tab or header) no corresponding
 * atom is given as an iterate; thus iterator may give 0, 1, or 2 iterates
 * depending on presence of adjacent atoms and specification of iteration
 * direction.
 */

public class AtomIteratorArrayListAdjacent implements AtomIteratorAtomDependent, java.io.Serializable {

    /**
     * Constructor gives iterator not ready for iteration.  Must
     * set an atom and reset before using.  Default direction is
     * null, meaning that both adjacent atoms (if there are two)
     * will be given by iterator.
     */
    public AtomIteratorArrayListAdjacent(IteratorDirective.Direction direction) {
        this(direction,new AtomToIndexChild(), new AtomToParentChildList());
    }
    
    public AtomIteratorArrayListAdjacent(IteratorDirective.Direction direction,
            AtomToIndex atomToIndex, AtomToArrayList atomToArrayList) {
        super();
        this.direction = direction;
        this.atomToIndex = atomToIndex;
        this.atomToArrayList = atomToArrayList;
        unset();
    }

    /**
     * Returns true if the given AtomSet has count == 1 and 
     * its atom is one of the iterates for the current condition
     * of the iterator (independent of hasNext status).
     */
    public boolean contains(AtomSet atom) {
        if(atom == null || atom.count() != 1) return false;
        Atom testAtom = atom.getAtom(0);
        if(testAtom == null) return false;
        if(direction != IteratorDirective.Direction.DOWN && firstCursor < list.size()-1 && list.get(firstCursor+1) == testAtom) {
            return true;
        }
        if(direction != IteratorDirective.Direction.UP && firstCursor > 0 && list.get(firstCursor-1) == testAtom) {
            return true;
        }
        return false;
    }

    /**
     * Returns the number of iterates that iterator would give
     * if reset and iterated in its current condition.  Does not
     * depend on or affect iteration state.
     */
    public int size() {
        int size = 0;
        if(direction != IteratorDirective.Direction.DOWN && (firstCursor < list.size()-1)) {
            size++;
        }
        if(direction != IteratorDirective.Direction.UP && (firstCursor > 0)) {
            size++;
        }
        return size;
    }

    /**
     * Performs action on all iterates for current condition of iterator.
     */
    public void allAtoms(AtomsetAction action) {
        if(direction != IteratorDirective.Direction.DOWN) {
            if (firstCursor < list.size()-1) {
                action.actionPerformed(list.get(firstCursor+1));
            }
        }
        if (direction != IteratorDirective.Direction.UP) {
            if (firstCursor > 0) {
                action.actionPerformed(list.get(firstCursor-1));
            }
        }
    }

    /**
     * Returns true if another iterate is forthcoming, false otherwise.
     */
    public boolean hasNext() {
        return hasNext;
    }

    /**
     * Returns 1, indicating that this is an atom AtomSet iterator.
     */
    public int nBody() {
        return 1;
    }

    /**
     * Same as nextAtom.
     */
    public AtomSet next() {
        return nextAtom();
    }

    /**
     * Returns the next iterator, or null if hasNext is false.
     */
    public Atom nextAtom() {
        if (!hasNext) {
            return null;
        }
        if(upListNow) {
            upListNow = false;
            if (direction == IteratorDirective.Direction.UP || firstCursor == 0) {
                hasNext = false;
            }
            return list.get(firstCursor+1);
        }
        hasNext = false;
        return list.get(firstCursor-1);
    }

    /**
     * Returns the next iterate without advancing the iterator.
     */
    public AtomSet peek() {
        if (!hasNext) {
            return null;
        }
        if(upListNow) {
            return list.get(firstCursor+1);
        }
        return list.get(firstCursor-1);
    }

    /**
     * Readies the iterator to begin iteration.
     */
    public void reset() {
        upListNow = false;
        hasNext = true;
        if ((direction != IteratorDirective.Direction.DOWN) && (firstCursor < list.size()-1)) {
            upListNow = true;
            return;
        }
        if ((direction != IteratorDirective.Direction.UP) && (firstCursor > 0)) {
            return;
        }
        hasNext = false;
    }

    /**
     * Puts iterator in a state where hasNext is false.
     */
    public void unset() {
        hasNext = false;
    }   
    
    /**
     * Sets the first atom for iteration. Iteration proceeds from this atom up
     * and/or down the list, depending on how iterator was configured at
     * construction.
     */
    public void setAtom(Atom atom) {
        list = atomToArrayList.getArrayList(atom);
        firstCursor = atomToIndex.getIndex(atom);
    }

    private int firstCursor;
    private final Direction direction;
    private boolean upListNow;
    private AtomArrayList list;
    private boolean hasNext;
    private final AtomToIndex atomToIndex;
    private final AtomToArrayList atomToArrayList;
}
