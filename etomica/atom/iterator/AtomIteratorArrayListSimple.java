package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;

 /**
  * An atom iterator of the elements from an AtomArrayList (in proper
  * sequence).  Iterator will fail if element are added to or removed 
  * from list while iteration is proceeding.
  */

public class AtomIteratorArrayListSimple implements AtomIterator, java.io.Serializable {

    /**
     * Constructs new iterator with an empty list.
     */
 	public AtomIteratorArrayListSimple() {
 		this(new AtomArrayList());
 	}
    
    /**
     * Constructs new iterator set to iterate given list (upon reset).
     */
 	public AtomIteratorArrayListSimple(AtomArrayList atomList) {
 		list = atomList;
 	}
    
    /**
     * Sets the list for iteration.  Null value is permitted, which will
     * cause iterator to give not iterates.
     */
 	public void setList(AtomArrayList atomList) {
        if(atomList != null) {
            list = atomList;
        } else {
            emptyList.clear();
            list = emptyList;
        }
 	}
 	
    /**
     * Returns 1, indicating that this is an atom iterator.
     */
 	public int nBody() {
        return 1;
    }
    
    /**
     * Puts iterator in state in which hasNext is false.
     */
 	public void unset() {
        cursor = size();
    }
 
    /**
     * Indicates if iterator has another iterate.
     */
 	public boolean hasNext() {
 	    return cursor < size();
 	}
 
    /**
     * Returns the next iterate and advances the iterator.
     */
 	public Atom nextAtom() {
        if (cursor < size()) {
            return list.get(cursor++);
        }
        return null;
 	}
 	
    /**
     * Same as nextAtom().
     */
 	public AtomSet next() {
 		return nextAtom();
 	}
 
    /**
     * Returns the next iterate without advancing the iterator.
     */
 	public AtomSet peek() {
 		return list.get(cursor);
 	}
    
    /**
     * Returns the number of iterates that would be given by this iterator
     * if reset with the current list.
     */
 	public int size() {
 		return list.size();
 	}

    /**
     * Performs action on all elements of current list.
     */
 	public void allAtoms(AtomsetAction act) {
 		int arraySize = size();
 		for (int i=0; i<arraySize; i++) {
 			act.actionPerformed(list.get(i));
 		}
 	}

    /**
     * Puts iterator in state ready to begin iteration.
     */
 	public void reset() {
 		cursor = 0;
 	}
 	
    /**
     * Returns true if the given atom is in the list.
     */
 	public boolean contains(AtomSet atom) {
        if(atom == null || atom.count() != 1) {
            return false;
        }
 		return list.contains(atom.getAtom(0));
 	}
 
    /**
     * Index of element to be returned by subsequent call to next.
     */
    protected int cursor = 0;

    protected AtomArrayList list;
    private final AtomArrayList emptyList = new AtomArrayList();
 }
