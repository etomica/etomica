package etomica;

import etomica.nbr.AtomSequencerNbr;

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
	
    /**
     * index of phase of interest.  -1 indicates no particular phase.
     */
    public static final int PHASE_INDEX = 0;
    
	/**
	 * true if debugging is currently enabled (when the integrator reaches step START) 
	 */
	public static boolean DEBUG_NOW = false;

	/**
	 * determines whether any of the atoms in the given array are set to be debugged
	 * (via ATOMx_NUM and setAtoms(phase)).
	 * @param atoms array of atoms to be checked for debugging status
	 * @return true if any of the atoms in the atoms array should be debugged
	 */
	public static boolean anyAtom(Atom[] atoms) {
		for (int i=0; i<atoms.length; i++) {
			if (atoms[i] == ATOM1 || atoms[i] == ATOM2) return true;
		}
		return false;
	}

	/**
	 * determines if all of the atoms in the given array are set to be debugged
	 * (via ATOMx_NUM and setAtoms(phase)).
	 * @param atoms array of atoms to be checked for debugging status
	 * @return true if all of the atoms in the atoms array should be debugged
	 */
	public static boolean allAtoms(Atom[] atoms) {
		for (int i=0; i<atoms.length; i++) {
			if (atoms[i] != ATOM1 && atoms[i] != ATOM2) return false;
		}
		return true;
	}

    /**
     * Checks whether the given phase is of debugging interest
     * @param checkPhase phase to be checked
     * @return true if the phase is of interest
     */
    public static boolean thisPhase(Phase checkPhase) {
         return checkPhase.index == PHASE_INDEX;
    }
    
	/**
	 * Atoms to be debugged.  These are set by setAtoms(phase)
	 */
	public static Atom ATOM1, ATOM2;

	/**
	 * Sets atoms to be debugged.  This sets ATOM1 and ATOM2 to be the
	 * ATOM1_NUMth and ATOM2_NUMth atoms in the phase
	 * @param phase the phase containing atoms to be debugged
	 */
	public static void setAtoms(Phase phase) {
		if (ATOM1_NUM > -1) ATOM1 = phase.speciesMaster.atomList.get(ATOM1_NUM);
		if (ATOM2_NUM > -1) ATOM2 = phase.speciesMaster.atomList.get(ATOM2_NUM);
	}

	/**
	 * Gives atoms ATOM1 and/or ATOM2 a thorough checking for sanity and consistency. 
	 * Computational cost and level of output depend on LEVEL. 
	 * @param cPair coordinate to use for calculating distance between atoms.
	 */
	public static void checkAtoms(Space.CoordinatePair cPair) {
		cPair.reset(Debug.ATOM1.coord, Debug.ATOM2.coord);
		double r2 = cPair.r2();
		//XXX What's the pair hard-core diameter?  Elephino!  Could check energy instead.  fun.
		if (Debug.LEVEL > 1 || Math.sqrt(r2) < Default.ATOM_SIZE-1.e-11) {
			System.out.println("distance between "+Debug.ATOM1+" and "+Debug.ATOM2+" is "+Math.sqrt(r2));
    		if (Debug.LEVEL > 2 || Math.sqrt(r2) < Default.ATOM_SIZE-1.e-11) {
    			System.out.println(Debug.ATOM1+" coordinates "+Debug.ATOM1.coord.position());
    			System.out.println(Debug.ATOM2+" coordinates "+Debug.ATOM2.coord.position());
    		}
/*			if (Math.sqrt(r2) < Default.ATOM_SIZE-1.e-11) {
				throw new RuntimeException("overlap");
			}*/
		}
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