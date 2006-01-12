package etomica.util;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomSet;
import etomica.phase.Phase;

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
	public static int START = 0;
	
	/**
	 * what step of the integrator should the simulation bail
	 * set this to prevent Etomica from running for a long time while 
	 * outputting debugging information (potentially filling up the disk)  
	 */
	public static int STOP = -1;
	
	/**
	 * debugging level.  A higher level means more debugging
	 * performance should suffer as the level goes up, and
	 * finding useful information might get difficult.
	 */
	public static int LEVEL = 1;
	
    /**
     * true if debugging is currently enabled (when the integrator reaches step START) 
     */
    public static boolean DEBUG_NOW = false;

	/**
	 * leaf atom number of first atom of interest.  The atom is the nth leaf atom 
	 * in a phase (set by calling setAtoms(phase)) More debugging information will be
	 * printed out about this particular atom.  -1 indicates no particular atom.
	 */
	public static final int ATOM1_NUM = -1;
	
	/**
	 * leaf atom number of second atom of interest.  This is often used in conjunction with 
	 * ATOM1_INDEX to collect information about a pair of atoms.  -1 indicates no
	 * particular atom.  
	 */
	public static final int ATOM2_NUM = -1;
	
    public static final int MOLECULE1_INDEX = -1;
    public static final int MOLECULE2_INDEX = -1;
    
    /**
     * index of phase of interest.  -1 indicates no particular phase.
     */
    public static final int PHASE_INDEX = 1;
    
    public static final double ATOM_SIZE = 1.0;
    
    /**
     * determines whether this atom is set to be debugged
     * (via ATOMx_NUM and setAtoms(phase)).
     * @param atom Atom to be checked for debugging status
     * @return true if the atom should be debugged
     */
    public static boolean thisAtom(Atom atom) {
        if (atom == ATOM1 || atom == ATOM2) return true;
        return false;
    }

	/**
	 * determines whether any of the atoms in the given array are set to be debugged
	 * (via ATOMx_NUM and setAtoms(phase)).
	 * @param atoms array of atoms to be checked for debugging status
	 * @return true if any of the atoms in the atoms array should be debugged
	 */
	public static boolean anyAtom(AtomSet atoms) {
		for (int i=0; i<atoms.count(); i++) {
			if ((ATOM1 != null && atoms.getAtom(i) == ATOM1) || (ATOM2 != null && atoms.getAtom(i) == ATOM2)) return true;
            if ((MOLECULE1 != null && atoms.getAtom(i) == MOLECULE1) || (MOLECULE2 != null && atoms.getAtom(i) == MOLECULE2)) return true;
		}
		return false;
	}

	/**
	 * determines if all of the atoms in the given array are set to be debugged
	 * (via ATOMx_NUM and setAtoms(phase)).
	 * @param atoms array of atoms to be checked for debugging status
	 * @return true if all of the atoms in the atoms array should be debugged
	 */
	public static boolean allAtoms(AtomSet atoms) {
		for (int i=0; i<atoms.count(); i++) {
			if (atoms.getAtom(i) != ATOM1 && atoms.getAtom(i) != ATOM2 
               && atoms.getAtom(i) != MOLECULE1 && atoms.getAtom(i) != MOLECULE2) return false;
		}
		return true;
	}

    /**
     * Checks whether the given phase is of debugging interest
     * @param checkPhase phase to be checked
     * @return true if the phase is of interest
     */
    public static boolean thisPhase(Phase checkPhase) {
         return checkPhase.getOrdinal() == PHASE_INDEX;
    }
    
	/**
	 * Atoms to be debugged.  These are set by setAtoms(phase)
	 */
	public static AtomLeaf ATOM1, ATOM2;
    public static Atom MOLECULE1, MOLECULE2;

	/**
	 * Sets atoms to be debugged.  This sets ATOM1 and ATOM2 to be the
	 * ATOM1_NUMth and ATOM2_NUMth atoms in the phase
	 * @param phase the phase containing atoms to be debugged
	 */
	public static void setAtoms(Phase phase) {
		if (ATOM1_NUM > -1) ATOM1 = (AtomLeaf)phase.getSpeciesMaster().leafList.get(ATOM1_NUM);
		if (ATOM2_NUM > -1) ATOM2 = (AtomLeaf)phase.getSpeciesMaster().leafList.get(ATOM2_NUM);
        if (MOLECULE1_INDEX > -1) MOLECULE1 = phase.molecule(MOLECULE1_INDEX);
        if (MOLECULE2_INDEX > -1) MOLECULE2 = phase.molecule(MOLECULE2_INDEX);
	}

}
