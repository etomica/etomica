/*
 * History
 * Created on Aug 23, 2004 by kofke
 */
package etomica.action;

import etomica.Atom;
import etomica.AtomSet;

/**
 * Action that simply detects if a particular set of atoms is ever
 * given as the argument of the actionPerformed method.  This is
 * considered to occur if a given array has elements that match,
 * in order, the elements of the array given to the constructor
 * at instantiation.
 */

 //used by some iterators to implement their contains() method

public class AtomsetDetect extends AtomActionAdapter {

	/**
	 * Constructs to test for the given atoms. Given array is copied, so
	 * test is always performed against the set of atoms in the array
	 * at time of instantiation. 
	 * @param testAtom atoms against which those passed to actionPerformed are compared.
	 */
	public AtomsetDetect(AtomSet testAtom) {
		atoms = testAtom;
	}
	
	/**
	 * Sets detect to true if atoms in the given array are the same 
	 * (in order) as the atoms specified via setAtoms.  Detect will stay true until reset(), 
	 * regardless of outcome of subsequent calls to this method.
	 */
	public void actionPerformed(AtomSet atom) {
		if(detected) return;
		detected = atom.equals(atoms);
	}
    
    /**
     * Same as actionPerformed(AtomSet), taking a single Atom.
     */
    public void actionPerformed(Atom atom) {
        actionPerformed((AtomSet)atom);
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
}
