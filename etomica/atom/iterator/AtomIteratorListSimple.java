package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.atom.AtomList;
import etomica.atom.AtomListTabbed;
import etomica.atom.AtomSet;

/**
 * Iterates through the elements of an untabbed AtomList instance. Iterates
 * uplist only, from beginning to end of list. Iterator functions correctly in
 * situations where elements are removed from list after they are returned by
 * iterator (but only for iteration via hasNext/next, not allAtoms).
 * 
 * @author David Kofke
 */

public final class AtomIteratorListSimple extends AtomIteratorAdapter implements
        AtomsetIteratorListDependent {

    /**
     * Constructs iterator with an empty list for iteration.
     */
    public AtomIteratorListSimple() {
        this(new AtomList());
    }

    /**
     * Constructs iterator for iteration over the given list. Subsequent call to
     * reset() is needed before beginning iteration.
     * 
     * @param list
     */
    public AtomIteratorListSimple(AtomList list) {
        super(new AtomIteratorSequence(IteratorDirective.Direction.UP));
        setList(list);
    }

    /**
     * Sets the (untabbed) list containing the atoms that will be returned by this
     * iterator. Call to reset() is needed before beginning iteration. If
     * argument is null, an empty list is created as the iterator's list.
     * 
     * @throws IllegalArgumentException if newList is an instance of AtomListTabbed
     */
    public void setList(AtomList newList) {
        if(newList instanceof AtomListTabbed) throw new IllegalArgumentException("Cannot iterate tabbed list with this iterator");
        list = (newList != null) ? newList : new AtomList();
        unset();
    }

    /**
     * @return the list used for iteration.
     */
    public AtomList getList() {
        return list;
    }

    public void reset() {
        ((AtomIteratorSequence) iterator).setFirst(list.header.next);
        super.reset();
    }

    /**
     * Returns true if the given atom is in the list of iterates, false
     * otherwise.
     */
    public boolean contains(AtomSet atom) {
        if(atom == null) return false;
        return list.contains(atom.getAtom(0));
    }

    public void allAtoms(AtomsetAction action) {
        ((AtomIteratorSequence) iterator).setFirst(list.header.next);
        super.allAtoms(action);
    }

    /**
     * Returns the total number of iterates that can be returned by this
     * iterator, for its current list basis.
     */
    public int size() {
        return list.size();
    }

    private AtomList list;

}//end of AtomIteratorListSimple

