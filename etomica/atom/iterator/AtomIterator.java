/* History
 * Created on Aug 4, 2004
 */
package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.atom.Atom;
import etomica.atom.AtomSet;

/**
 * Interface for classes that loop over a set of atoms. Permits
 * iteration via a hasNext()-next() while loop (iterator returns
 * atoms to client) or via a call to allAtoms(AtomActive) (client gives
 * action to iterator).
 */

public interface AtomIterator extends AtomsetIterator {
                    
	/**
	 * Returns the next atom in the iteration sequence, or
     * null if hasNext() is false.
	 */
    public Atom nextAtom();
        
    /**
     * Static iterator that returns no atoms.
     * @author kofke
     */
    public static final AtomIterator NULL = new AtomIterator() {
    	public void allAtoms(AtomsetAction action) {}
    	public boolean contains(AtomSet atom) {return false;}
    	public boolean hasNext() {return false;}
       	public AtomSet next() {return null;}
       	public Atom nextAtom() {return null;}
    	public void reset() {}
    	public int size() {return 0;}
    	public AtomSet peek() {return null;}
    	public void unset() {}
    	public int nBody() {return 1;}
    };
}
