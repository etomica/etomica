/*
 * History
 * Created on Aug 23, 2004 by kofke
 */
package etomica;

/**
 * Action that simply detects if a particular atom is ever
 * given as the argument of the actionPerformed method.
 */

 //used by some iterators to implement their contains() method

public class AtomActiveDetect implements AtomActive {

	/**
	 * Constructs to test for the given atom.  null value is permitted,
	 * and will cause to test for nulls given to actionPerformed. 
	 * @param testAtom atom against which those passed to actionPerformed are compared.
	 * 
	 */
	public AtomActiveDetect(Atom testAtom) {
		this.testAtom = testAtom;
	}
	/**
	 * Sets detect to true if given atom is same as atom
	 * given in constructor.
	 */
	public void actionPerformed(Atom atom) {
		if(testAtom == null) detected |= (atom == null);
		else detected |= (testAtom.equals(atom)); //same as detected = (detected || testAtom.equals(atom))
	}
	
	/**
	 * Sets the callCount to zero.
	 */
	public void reset() {detected = false;}
	
	/**
	 * @return true if the testAtom was ever given as the
	 * argument to actionPerformed since the last call to
	 * reset (or instantation if reset was not yet called).
	 */
	public boolean detectedAtom() {return detected;}

	private boolean detected = false;
	private final Atom testAtom;
}
