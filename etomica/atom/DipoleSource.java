package etomica.atom;

import etomica.api.IMolecule;
import etomica.api.IVector;

/**
 * Interface for something that can calculate the dipole of a molecule.
 */
public interface DipoleSource {
    
    /**
     * Returns the dipole of the given molecule.
     * The method is likely to through an exception if a molecule of the wrong
     * type is passed in.
     */
    public IVector getDipole(IMolecule molecule);
}