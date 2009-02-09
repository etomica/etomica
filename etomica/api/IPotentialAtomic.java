package etomica.api;

public interface IPotentialAtomic extends IPotential {

    /**
     * Returns the interaction energy between the given atoms.  There might be
     * 0, 1, 2 or more atoms in the AtomSet.
     */
    public double energy(IAtomList atoms);
}
