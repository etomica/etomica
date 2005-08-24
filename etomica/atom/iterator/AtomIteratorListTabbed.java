package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.action.AtomsetDetect;
import etomica.atom.Atom;
import etomica.atom.AtomLinker;
import etomica.atom.AtomList;
import etomica.atom.AtomSet;

/**
 * Atom iterator that traverses the elements of a tabbed (or untabbed) atom
 * list. Configurable to permit iteration beginning at any specified tab, and
 * ending at any specified tab type. Iteration proceeds uplist only.
 * 
 * @see AtomListTabbed
 */

public final class AtomIteratorListTabbed implements AtomIterator,
        AtomsetIteratorListDependent, java.io.Serializable {

    /**
     * Constructs a new iterator using an empty list as its basis for iteration.
     */
    public AtomIteratorListTabbed() {
        this(new AtomList());
    }

    /**
     * Constructs a new iterator using the given list as its basis for
     * iteration. Iterator is conditioned to start from the first element of
     * list and iterate to the end, ignoring any tabs before reaching the end.
     * Iterator must be reset before use.
     */
    public AtomIteratorListTabbed(AtomList list) {
        sequenceIterator = new AtomIteratorSequence(IteratorDirective.UP);
        setList(list);
    }

    /**
     * Sets the given list of atoms as the basis for iteration. The atoms
     * returned by this iterator will be those from the given list. After
     * calling this method, iterator must be reset before use. Sets first and
     * terminator iteration elements to be the header of the given list,
     * regardless of the previous settings. If given a null argument, an empty
     * list is created for iteration.
     */
    public void setList(AtomList newList) {
        list = (newList != null) ? newList : new AtomList();
        startTab = list.header;
        terminatorType = AtomLinker.Tab.HEADER_TAB;
    }

    /**
     * @return the list defining the atoms given by this iterator.
     */
    public AtomList getList() {
        return list;
    }

    /**
     * Set the iterator to begin iteration in its current condition (first,
     * terminator type most recently set).
     */
    public void reset() {
        sequenceIterator.setFirst(startTab.next);
        sequenceIterator.reset();
        while (sequenceIterator.peekLinker().atom == null) {
            if (isTerminator((AtomLinker.Tab) sequenceIterator.nextLinker())) {
                unset();
                break;
            }
        }
    }

    /**
     * Same as nextAtom, for AtomsetIterator interface.
     */
    public final AtomSet next() {
        return nextAtom();
    }

    /**
     * Returns the next atom and advances the iterator. Returns null if hasNext
     * is false.
     */
    public Atom nextAtom() {
        Atom atom = sequenceIterator.nextAtom();
        while (sequenceIterator.peek() == null)
            if (isTerminator((AtomLinker.Tab) sequenceIterator.nextLinker())) {
                unset();
                break;
            }
        return atom;
    }

    /**
     * Returns the next iterate without advancing the iterator. Returns null if
     * hasNext is false.
     */
    public AtomSet peek() {
        return sequenceIterator.peek();
    }

    /**
     * Indicates whether iterator has another iterate.
     */
    public boolean hasNext() {
        return sequenceIterator.hasNext();
    }

    /**
     * Puts iterator in state in which hasNext is false.
     */
    public void unset() {
        sequenceIterator.unset();
    }

    /**
     * Returns 1, indicating that this is an atom AtomsetIterator
     */
    public int nBody() {
        return 1;
    }

    /**
     * Performs action on all atoms according to the current condition of the
     * iterator (first, terminator type). Does not require reset before use.
     * Clobbers iteration state.
     */
    public void allAtoms(AtomsetAction action) {
        sequenceIterator.setFirst(startTab.next);
        sequenceIterator.reset();
        do {
            AtomLinker link = sequenceIterator.nextLinker();
            if (link.atom != null) {
                action.actionPerformed(link.atom);
            } else if (isTerminator((AtomLinker.Tab) link)) {
                break;
            }
        } while (true);
    }//end of allAtoms

    //    /**
    //     * Sets iteration to be in the given direction and unsets iterator.
    //     * A null argument has the same effect as UP.
    //     */
    //    public void setDirection(IteratorDirective.Direction direction) {
    //        upList = (direction == IteratorDirective.UP) || (direction == null);
    //		unset();
    //    }
    //    
    //    /**
    //     * @return the current setting for the iteration direction.
    //     */
    //    public IteratorDirective.Direction getDirection() {
    //    	return upList ? IteratorDirective.UP : IteratorDirective.DOWN;
    //    }

    /**
     * Resets to begin with the given atom linker and unsets iterator. If given
     * linker is null, sets first to header of list.
     * 
     * @throws IllegalArgumentException
     *             if linker is not in the list currently set for iteration
     */
    public void setFirst(AtomLinker.Tab first) {
        if (first == null) {
            startTab = list.header;
        } else if (first.list != this.list) {
            throw new IllegalArgumentException(
                    "Attempt to iterate from a tab not in the list");
        } else {
            startTab = first;
        }
        unset();
    }

    private boolean isTerminator(AtomLinker.Tab tab) {
        return (tab.type & terminatorType) != 0;
    }

    /**
     * Sets the tab type of the terminator. Iteration stops when the header or a
     * tab of this type is encountered. A value of zero will cause only the
     * header to serve as the terminator; a negative value will cause any tab to
     * act as a terminator. Otherwise value is normally obtained from the type
     * field of a Tab instance.
     */
    public void setTerminatorType(int terminatorType) {
        if (terminatorType < 0) {
            this.terminatorType = AtomLinker.Tab.ANY_TAB;
        } else {
            this.terminatorType = terminatorType | AtomLinker.Tab.HEADER_TAB;
        }
    }

    /**
     * @return the current tab type that indicates termination of iteration
     */
    public int getTerminatorType() {
        return terminatorType;
    }

    /**
     * Returns true if the given atoms set is in the list of iterates that would
     * be given by iterator as currently conditioned; returns false otherwise.
     * Ignores any atoms other than the first in the array. Returns false if
     * array does not specify an atom, or if the count of the atomset is not
     * equal to 1. Does not require iterator be reset.
     */
    public boolean contains(AtomSet atom) {
        if (atom == null || atom.count() != 1) {
            return false;
        }
        AtomsetDetect detector = new AtomsetDetect(atom);
        allAtoms(detector);
        return detector.detectedAtom();
    }

    /**
     * Returns the total number of iterates that can be returned by this
     * iterator, for its current list basis, and as it is currently conditioned
     * (as given by most recently specified first and terminatorType). Does not
     * require the iterator be reset.
     */
    public int size() {
        if(startTab == list.header && terminatorType == AtomLinker.Tab.HEADER_TAB) {
            return list.size();
        }
        AtomsetCount counter = new AtomsetCount();
        allAtoms(counter);
        return counter.callCount();
    }

    private AtomList list;
    private int terminatorType;//bitmask to match against for terminators
    private final AtomIteratorSequence sequenceIterator;
    private AtomLinker.Tab startTab;

}//end of AtomIteratorList

