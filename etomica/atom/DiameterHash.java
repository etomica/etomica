package etomica.atom;

import etomica.api.IAtom;

/**
 * This provides an interface for a class that acts like a hashmap for atomic
 * diameters.
 *
 * @author Andrew Schultz
 */
public interface DiameterHash {

    /**
     * Returns an appropriate diameter for the given atom.  If no diameter is
     * known for the atom, -1 is returned.
     */
    public double getDiameter(IAtom atom);

}
