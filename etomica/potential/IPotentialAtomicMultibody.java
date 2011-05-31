package etomica.potential;

import etomica.api.IPotentialAtomic;

/**
 * Interface for an atomic non-additive potential.  Although the potential is
 * varies with the angles between the atoms, the potential is spherically
 * symmetric in that the atoms do not have orientation and the potential can be
 * computed from the iteratomic distances.
 *
 * @author Andrew Schultz
 */
public interface IPotentialAtomicMultibody extends IPotentialAtomic {
    public double energy(double[] r2);
}