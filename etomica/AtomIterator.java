/* History
 * Created on Aug 4, 2004
 */
package etomica;

import etomica.action.AtomsetAction;

/**
 * Interface for classes that loop over a set of atoms. Permits
 * iteration via a hasNext()-next() while loop (iterator returns
 * atoms to client) or via a call to allAtoms(AtomActive) (client gives
 * action to iterator).
 */

public interface AtomIterator extends AtomsetIterator {
                    
	/**
	 * Returns the next atom in the iteration sequence.
	 * No specific behavior is guaranteed if hasNext() == false 
	 * at the time method is called, except that calling next()
	 * will not cause hasNext() to become true.
	 */
    public Atom nextAtom();
        
    /**
     * Static iterator that returns no atoms.
     * @author kofke
     */
    public static final AtomIterator NULL = new AtomIterator() {
    	private final Atom[] atoms = new Atom[1];
    	public void allAtoms(AtomsetAction action) {}
    	public boolean contains(Atom[] atom) {return false;}
     	public boolean contains(Atom atom) {return false;}
    	public boolean hasNext() {return false;}
       	public Atom[] next() {atoms[0] = null; return atoms;}
       	public Atom nextAtom() {return null;}
    	public void reset() {}
    	public int size() {return 0;}
    	public Atom[] peek() {atoms[0] = null; return atoms;}
    	public void unset() {}
    	public int nBody() {return 1;}
    };
}
