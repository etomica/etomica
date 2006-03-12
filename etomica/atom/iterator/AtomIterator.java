/* History
 * Created on Aug 4, 2004
 */
package etomica.atom.iterator;

import etomica.atom.Atom;

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
}
