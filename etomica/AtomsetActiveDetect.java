/*
 * History
 * Created on Aug 23, 2004 by kofke
 */
package etomica;

/**
 * Action that simply detects if a particular set of atoms is ever
 * given as the argument of the actionPerformed method.  This is
 * consider to occur if a given array has elements that match,
 * in order, the elements of the array given to the constructor
 * at instantiation.
 */

 //used by some iterators to implement their contains() method

public class AtomsetActiveDetect implements AtomsetActive {

	/**
	 * Constructs to test for the given atoms. Given array is copied, so
	 * test is always performed against the set of atoms in the array
	 * at time of instantiation. 
	 * @param testAtom atoms against which those passed to actionPerformed are compared.
	 */
	public AtomsetActiveDetect(Atom[] testAtom) {
		this.testAtom = (Atom[])testAtom.clone();
	}
	
	/**
	 * Constructs to test for the given atom.
	 * @param testAtom atom against which those passed to actionPerformed are compared.
	 */
	public AtomsetActiveDetect(Atom testAtom) {
		this.testAtom = new Atom[] {testAtom};
	}
	
	/**
	 * Sets detect to true if atoms in the given array are the same 
	 * (in order) as the atoms given in constructor at time of 
	 * instantiation.  Detect will stay true until reset(), 
	 * regardless of outcome of subsequent calls to this method.
	 */
	public void actionPerformed(Atom[] atom) {
		if(detected) return;
		detected = java.util.Arrays.equals(atom, testAtom);
	}
	
	/**
	 * Sets detected to false.
	 */
	public void reset() {detected = false;}
	
	/**
	 * @return true if the testAtom was ever given as the
	 * argument to actionPerformed since the last call to
	 * reset (or instantation if reset was not yet called).
	 */
	public boolean detectedAtom() {return detected;}

	private boolean detected = false;
	private final Atom[] testAtom;
}
