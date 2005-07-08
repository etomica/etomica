package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomIterator;
import etomica.AtomSet;
import etomica.action.AtomsetAction;
import etomica.atom.AtomList;
import etomica.utility.Arrays;

/**
 * Iterates over all the atoms given across an array of atom lists. Atoms
 * present in multiple lists, or in lists added multiple times to iterator, will
 * be iterated multiple times.
 */
public class AtomIteratorListCompound implements AtomIterator, java.io.Serializable {

    /**
     * Construct iterator to loop over all atoms in each list in the given
     * array. Must reset before using.  Given array is cloned for internal
     * use, but references to lists it contains are maintained (lists are
     * not copied).
     */
    public AtomIteratorListCompound(AtomList[] lists) {
        listSet = (AtomList[])lists.clone();
        iterator = new AtomIteratorListSimple();
        unset();
    }

    /**
     * Constructs iterator that will give no iterates until a list is added.
     */
    public AtomIteratorListCompound() {
        this(new AtomList[0]);
    }

    /**
     * Adds the given list to the set of lists collected by this iterator. No
     * check is made that the list was not already added; if a list is added
     * more than once, it will be iterated for each addition, as if it were a
     * different instance each time. If list is null, no action is taken. 
     * Iterator is left in an unset state.
     */
    public void addList(AtomList list) {
        if (list == null) {
            return;
        }
        listSet = (AtomList[]) Arrays.addObject(listSet, list);
        unset();
    }

    /**
     * Sets the lists of atoms that will be iterated, discarding any lists
     * previously added. If argument is null, no iterates will be given.
     * Reference to lists is maintained, so changes to them will be reflected by
     * the iterator; reference to array is not kept (clone of array is used
     * internally by iterator). Leaves iterator unset.
     */
    public void setLists(AtomList[] lists) {
        if (lists == null) {
            listSet = new AtomList[0];
        } else {
            listSet = (AtomList[]) lists.clone();
        }
        unset();
    }

    /**
     * Removes the given list from the set of lists collected by this iterator.
     * If list was not previously added (or is null), no action is taken. If
     * list was added multiple times, removes only one reference to it. Leaves
     * iterator unset.
     */
    public void removeList(AtomList list) {
        listSet = (AtomList[]) Arrays.removeObject(listSet, list);
        unset();
    }

    /**
     * Indicates whether iterator has another iterate.
     */
    public boolean hasNext() {
        return iterator.hasNext();
    }

    /**
     * Puts iterator in state in which hasNext is false.
     */
    public void unset() {
        iterator.unset();
    }

    /**
     * Returns 1, indicating that this is an atom iterator.
     */
    public int nBody() {
        return 1;
    }

    /**
     * Returns the number of iterates that would be given on reset of iterator.
     * Equal to the sum of the sizes of all lists added to iterator.
     */
    public int size() {
        int count = 0;
        for (int i = 0; i < listSet.length; i++) {
            count += listSet[i].size();
        }
        return count;
    }

    /**
     * Returns true if the given atom is among the iterates given by this
     * iterator. Does not affect iteration state. Returns false if argument is
     * null or if atoms.count != 1.
     */
    public boolean contains(AtomSet atoms) {
        if (atoms == null || atoms.count() != 1) {
            return false;
        }
        for (int i = 0; i < listSet.length; i++) {
            if (listSet[i].contains(atoms.getAtom(0)))
                return true;
        }
        return false;
    }

    /**
     * Readies iterator to begin iteration.
     */
    public void reset() {
        if (listSet.length == 0) {
            return;
        }
        index = 0;
        while (listSet[index].isEmpty() && index + 1 < listSet.length) {
            index++;
        }
        if (index < listSet.length) {
            iterator.setList(listSet[index]);
            iterator.reset();
        }
    }

    /**
     * Returns next iterate without advancing the iterator.
     */
    public AtomSet peek() {
        return iterator.peek();
    }

    /**
     * Same as nextAtom.
     */
    public AtomSet next() {
        return nextAtom();
    }

    /**
     * Returns next iterate and advances the iterator.
     */
    public Atom nextAtom() {
        if (!hasNext()) {
            return null;
        }
        Atom atom = iterator.nextAtom();
        while (!iterator.hasNext()) {
            if (++index < listSet.length) {
                iterator.setList(listSet[index]);
                iterator.reset();
            } else {
                break;
            }
        }
        return atom;
    }

    /**
     * Performs action on all atoms in added lists. Performs action multiple
     * times to atoms present in multiple lists, or to lists added multiple
     * times to iterator.
     */
    public void allAtoms(AtomsetAction action) {
        for (int i = 0; i < listSet.length; i++) {
            iterator.setList(listSet[i]);
            iterator.reset();
            iterator.allAtoms(action);
        }
    }

    private AtomList[] listSet;
    private final AtomIteratorListSimple iterator;
    private int index;

}//end of AtomIteratorListCompound
