package etomica.util;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.integrator.IntegratorHard;
import etomica.nbr.list.AtomSequencerNbr;
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
         return checkPhase.getIndex() == PHASE_INDEX;
    }
    
	/**
	 * Atoms to be debugged.  These are set by setAtoms(phase)
	 */
	public static Atom ATOM1, ATOM2;
    public static Atom MOLECULE1, MOLECULE2;

	/**
	 * Sets atoms to be debugged.  This sets ATOM1 and ATOM2 to be the
	 * ATOM1_NUMth and ATOM2_NUMth atoms in the phase
	 * @param phase the phase containing atoms to be debugged
	 */
	public static void setAtoms(Phase phase) {
		if (ATOM1_NUM > -1) ATOM1 = phase.getSpeciesMaster().atomList.get(ATOM1_NUM);
		if (ATOM2_NUM > -1) ATOM2 = phase.getSpeciesMaster().atomList.get(ATOM2_NUM);
        if (MOLECULE1_INDEX > -1) MOLECULE1 = phase.molecule(MOLECULE1_INDEX);
        if (MOLECULE2_INDEX > -1) MOLECULE2 = phase.molecule(MOLECULE2_INDEX);
	}

	/**
	 * Gives atoms ATOM1 and/or ATOM2 a thorough checking for sanity and consistency. 
	 * Computational cost and level of output depend on LEVEL. 
	 * @param cPair coordinate to use for calculating distance between atoms.
	 */
	public static void checkAtoms() {
		if (Debug.ATOM1.seq instanceof AtomSequencerNbr) {
			int dnNbrCount = 0, upNbrCount = 0;
			AtomArrayList[] list = ((AtomSequencerNbr)Debug.ATOM1.seq).getDownList();
			for (int i=0; i<list.length; i++) {
				if (list[i].contains(Debug.ATOM2)) {
					if (Debug.LEVEL > 2) System.out.println(Debug.ATOM2+" is a down neighbor of "+Debug.ATOM1);
					dnNbrCount++;
				}
			}
			if (dnNbrCount > 1) throw new RuntimeException(Debug.ATOM2+" is a neighbor of "+Debug.ATOM1+" more than once");
			list = ((AtomSequencerNbr)Debug.ATOM1.seq).getUpList();
			for (int i=0; i<list.length; i++) {
				if (list[i].contains(Debug.ATOM2)) {
					if (Debug.LEVEL > 2) System.out.println(Debug.ATOM2+" is an up neighbor of "+Debug.ATOM1);
       				if (dnNbrCount > 0) throw new RuntimeException(Debug.ATOM2+" is a down and up neighbor of "+Debug.ATOM1);
					upNbrCount++;
				}
			}
            if (ATOM1.ia instanceof IntegratorHard.Agent && upNbrCount == 0) {
                if (((IntegratorHard.Agent)ATOM1.ia).collisionPartner == ATOM2) {
                    throw new IllegalStateException(ATOM2+" is collision partner of "+ATOM1+" but isn't in "+ATOM1+"'s uplist of neighbors");
                }
            }
			list = ((AtomSequencerNbr)Debug.ATOM2.seq).getDownList();
			for (int i=0; i<list.length; i++) {
				if (list[i].contains(Debug.ATOM1)) {
					if (Debug.LEVEL > 2) System.out.println(Debug.ATOM1+" is a down neighbor of "+Debug.ATOM2);
       				if (dnNbrCount > 0) throw new RuntimeException("too much down neighborness between "+Debug.ATOM2+" and "+Debug.ATOM1);
       				if (upNbrCount == 0) {
       					throw new RuntimeException(Debug.ATOM1+" is a down neighbor of "+Debug.ATOM2+" but not the other way around");
       				}
					dnNbrCount++;
				}
			}
			list = ((AtomSequencerNbr)Debug.ATOM2.seq).getUpList();
			for (int i=0; i<list.length; i++) {
				if (list[i].contains(Debug.ATOM1)) {
					if (Debug.LEVEL > 2) System.out.println(Debug.ATOM2+" is an up neighbor of "+Debug.ATOM1);
       				if (upNbrCount > 0) throw new RuntimeException("too much up neighborness between "+Debug.ATOM2+" and "+Debug.ATOM1);
       				if (dnNbrCount == 0) {
       					throw new RuntimeException(Debug.ATOM2+" is an up neighbor of "+Debug.ATOM1+" but not the other way around");
       				}
					upNbrCount++;
				}
			}
		}
	}

}
