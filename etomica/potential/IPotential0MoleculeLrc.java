package etomica.potential;

import etomica.api.IMolecule;
import etomica.api.IPotential;

/**
 * Interface for a long-range correction potential.
 *
 * @author Andrew Schultz
 */
public interface IPotential0MoleculeLrc extends IPotential {

    /**
     * Informs the potential of a target atom.  Null target atom indicates no
     * target atom.
     */
    public void setTargetMolecule(IMolecule targetAtom);

}