/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSetSinglet;

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
 	public AtomIteratorArrayListSimple(IAtomList atomList) {
 		list = atomList;
        atomSetSinglet = new AtomSetSinglet();
 	}
    
    /**
     * Sets the list for iteration.  Null value will result in a
     * NullPointerException.
     */
 	public void setList(IAtomList atomList) {
        list = atomList;
 	}
 	
 	public IAtomList getList() {
        return list;
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
        cursor = list.getAtomCount();
    }
 
    /**
     * Returns the next iterate and advances the iterator.
     */
 	public IAtom nextAtom() {
        if (cursor < list.getAtomCount()) {
            return list.getAtom(cursor++);
        }
        return null;
 	}
 	
    /**
     * Same as nextAtom().
     */
 	public IAtomList next() {
        IAtom atom = nextAtom();
        if (atom == null) return null;
        atomSetSinglet.atom = atom;
 		return atomSetSinglet;
 	}
 
    /**
     * Returns the number of iterates that would be given by this iterator
     * if reset with the current list.
     */
 	public int size() {
 		return list.getAtomCount();
 	}

    /**
     * Puts iterator in state ready to begin iteration.
     */
 	public void reset() {
 		cursor = 0;
 	}
 	
    private static final long serialVersionUID = 2L;

    /**
     * Index of element to be returned by subsequent call to next.
     */
    protected int cursor = 0;

    protected IAtomList list;
    protected final AtomSetSinglet atomSetSinglet;
 }
