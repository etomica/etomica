package etomica;

/**
 * Class holding static fields that determine whether debugging is on, how
 * much, what types and what (if anything) should be looked at specifically.
 * @author andrew
 */

public final class Debug {

	/**
	 * true if any debugging should be done
	 */
	public static final boolean ON = false;
        
	/**
	 * what step of the integrator deubgging should start from.
	 * to get debugging to start before Integrator.run, explicitly
	 * set DEBUG_NOW to true. 
	 */
	public static final int START = 0;
	
	/**
	 * what step of the integrator should the simulation bail
	 * set this to prevent Etomica from running for a long time while 
	 * outputting debugging information (potentially filling up the disk)  
	 */
	public static final int STOP = -1;
	
	/**
	 * debugging level.  A higher level means more debugging
	 * performance should suffer as the level goes up, and
	 * finding useful information might get difficult.
	 */
	public static final int LEVEL = 1;
	
	/**
	 * index of first atom of interest.  More debugging information will be
	 * printed out about this particular atom.  -1 indicates no particular atom.
	 */
	public static final int ATOM1_INDEX = -1;
	
	/**
	 * index of second atom of interest.  This is often used in conjunction with 
	 * ATOM1_INDEX to collect information about a pair of atoms.  -1 indicates no
	 * particular atom.  
	 */
	public static final int ATOM2_INDEX = -1;
	
	/**
	 * true if debuggin is currently enabled (when the integrator reaches step START) 
	 */
	public static boolean DEBUG_NOW = true;

	/**
	 * determines whether any of the atoms in the given array are set to be debugged
	 * (via ATOMx_INDEX). 
	 * @param atoms array of atoms to be checked for debugging status
	 * @return true if any of the atoms in the atoms array should be debugged
	 */
	public static boolean anyAtom(Atom[] atoms) {
		for (int i=0; i<atoms.length; i++) {
			if (atoms[i].node.index() == ATOM1_INDEX || atoms[i].node.index() == ATOM2_INDEX) return true;
		}
		return false;
	}

	/**
	 * determines if all of the atoms in the given array are set to be debugged
	 * (via ATOMx_INDEX). 
	 * @param atoms array of atoms to be checked for debugging status
	 * @return true if all of the atoms in the atoms array should be debugged
	 */
	public static boolean allAtoms(Atom[] atoms) {
		for (int i=0; i<atoms.length; i++) {
			if (atoms[i].node.index() != ATOM1_INDEX && atoms[i].node.index() != ATOM2_INDEX) return false;
		}
		return true;
	}

}