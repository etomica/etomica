package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.atom.AtomToAtomSet;
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
            AtomToIndex atomToIndex, AtomToAtomSet atomToAtomSet) {
        if (direction == null)
            throw new IllegalArgumentException(
                    "Must specify direction to constructor of AtomLinkerIterator");
        upListNow = (direction == IteratorDirective.Direction.UP);

        this.atomToIndex = atomToIndex;
        this.atomToAtomSet = atomToAtomSet;
        this.numToSkip = numToSkip;
 	}
    
    /**
     * Returns the next iterate and advances the iterator.
     */
 	public IAtom nextAtom() {
        if (upListNow) {
            if (cursor > list.getAtomCount()-2) {
                return null;
            }
            cursor++;
        }
        else {
            if (cursor < 1) {
                return null;
            }
            cursor--;
        }
 		return list.getAtom(cursor);
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
        IAtomSet localList = atomToAtomSet.getAtomSet(startAtom);
        int firstCursor = atomToIndex.getIndex(startAtom);
 		int arraySize = localList.getAtomCount();
        if (upListNow) {
            for (int i=firstCursor+numToSkip; i<arraySize; i++) {
                atomSetSinglet.atom = localList.getAtom(i);
                act.actionPerformed(atomSetSinglet);
            }
        }
        else {
            for (int i=firstCursor-numToSkip; i>-1; i--) {
                atomSetSinglet.atom = localList.getAtom(i);
                act.actionPerformed(atomSetSinglet);
            }
        }
 	}

    /**
     * Puts iterator in state ready to begin iteration.
     */
 	public void reset() {
        list = atomToAtomSet.getAtomSet(startAtom);
 		cursor = atomToIndex.getIndex(startAtom);
        if (upListNow) {
            cursor--;
            cursor += numToSkip;
        }
        else {
            cursor++;
            cursor -= numToSkip;
        }
 	}
    
    public void unset() {
        if (upListNow) {
            cursor = list.getAtomCount();
        }
        else {
            cursor = -1;
        }
    }
 	
    /**
     * Sets the first atom for iteration. Iteration proceeds from this atom up
     * and/or down the list, depending on how iterator was configured at
     * construction.
     */
    public void setAtom(IAtom atom) {
        startAtom = atom;
    }

    private static final long serialVersionUID = 1L;
    protected final boolean upListNow;
    private final AtomsetCount counter = new AtomsetCount();
    private final int numToSkip;
    private IAtom startAtom;
    private final AtomToIndex atomToIndex;
    private final AtomToAtomSet atomToAtomSet;
 }
