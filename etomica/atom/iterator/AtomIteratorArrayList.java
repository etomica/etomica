package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.action.AtomsetDetect;
import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.atom.AtomToArrayList;
import etomica.atom.AtomToIndex;
import etomica.atom.AtomToIndexChild;
import etomica.atom.AtomToParentChildList;

 /**
  * An atom iterator of the elements from an AtomArrayList (in proper
  * sequence).  Iterator will fail if element are added to or removed 
  * from list while iteration is proceeding.
  */
public class AtomIteratorArrayList extends AtomIteratorArrayListSimple implements AtomIteratorAtomDependent {

    /**
     * Constructs new iterator with an empty list.
     */
 	public AtomIteratorArrayList(IteratorDirective.Direction direction) {
 		this(direction, 0);
 	}
    
    /**
     * Constructs new iterator set to iterate given list (upon reset).
     */
 	public AtomIteratorArrayList(IteratorDirective.Direction direction, int numToSkip) {
        this(direction, numToSkip, new AtomToIndexChild(), new AtomToParentChildList());
    }
    
    public AtomIteratorArrayList(IteratorDirective.Direction direction, int numToSkip, 
            AtomToIndex atomToIndex, AtomToArrayList atomToArrayList) {
        if (direction == null)
            throw new IllegalArgumentException(
                    "Must specify direction to constructor of AtomLinkerIterator");
        upListNow = (direction == IteratorDirective.Direction.UP);

        this.atomToIndex = atomToIndex;
        this.atomToArrayList = atomToArrayList;
        this.numToSkip = numToSkip;
 	}
    
    /**
     * Returns the next iterate and advances the iterator.
     */
 	public Atom nextAtom() {
        if (cursor == -1) {
            return null;
        }
        int oldCursor = cursor;
        if (upListNow) {
            cursor++;
            if (cursor == list.size()) {
                cursor = -1;
            }
        }
        else {
            cursor--;
        }
 		return list.get(oldCursor);
 	}
 	
    /**
     * Returns the number of iterates that would be given by this iterator
     * if reset with the current list.
     */
    public int size() {
        counter.reset();
        allAtoms(counter);
        return counter.callCount();
    }

    /**
     * Performs action on all elements of current list.
     */
 	public void allAtoms(AtomsetAction act) {
 		int arraySize = list.size();
        if (upListNow) {
            for (int i=firstCursor+numToSkip; i<arraySize; i++) {
                act.actionPerformed(list.get(i));
            }
        }
        else {
            for (int i=firstCursor-numToSkip; i>-1; i--) {
                act.actionPerformed(list.get(i));
            }
        }
 	}

    /**
     * Puts iterator in state ready to begin iteration.
     */
 	public void reset() {
 		cursor = firstCursor;
        if (upListNow) {
            cursor += numToSkip;
            if (cursor > list.size()-1) {
                cursor = -1;
            }
        }
        else {
            cursor -= numToSkip;
            if (cursor < 0) {
                cursor = -1;
            }
        }
 	}
    
    public boolean hasNext() {
        return cursor != -1;
    }
    
    public void unset() {
        cursor = -1;
    }
 	
    /**
     * Returns true if the given atom is in the list.
     */
 	public boolean contains(AtomSet atom) {
        if(atom == null || atom.count() != 1) {
            return false;
        }
        detector.setAtoms(atom);
        detector.reset();
        allAtoms(detector);
        return detector.detectedAtom();
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

    protected final boolean upListNow;
    private final AtomsetCount counter = new AtomsetCount();
    private final AtomsetDetect detector = new AtomsetDetect(null);
    private final int numToSkip;
    private int firstCursor;
    private final AtomToIndex atomToIndex;
    private final AtomToArrayList atomToArrayList;
 }
