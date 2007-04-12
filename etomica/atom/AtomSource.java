package etomica.atom;

import etomica.phase.Phase;

/**
 * Interface for objects when return atoms (meeting some specification)
 * from a phase.
 */
public interface AtomSource {
    
    /**
     * sets the Phase the source should pull Atoms from.
     * Phase should not be null
     */
    public void setPhase(Phase p);

    /**
     * Returns an atom.  Will return null if there are no appropriate atoms in 
     * the given phase.
     */
    public IAtom getAtom();
}
