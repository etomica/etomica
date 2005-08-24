/* History
 * Created on Aug 4, 2004
 */
package etomica.atom;


/**
 * Interface for classes that loop over pairs of atoms. Permits
 * iteration via a hasNext()-nextPair() while loop.
 */

public interface AtomPairIterator extends AtomsetIterator {
                    
	/**
	 * Returns the next atom in the iteration sequence.
	 * No specific behavior is guaranteed if hasNext() == false 
	 * at the time method is called, except that calling next()
	 * will not cause hasNext() to become true.
	 */
    public AtomPair nextPair();
        
}
