package etomica.potential;

import etomica.api.IAtomType;

/**
 * Wraps a soft-spherical potential to apply a truncation to it.  Energy and
 * its derivatives are set to zero at a specified cutoff.  (No accounting is
 * made of the infinite force existing at the cutoff point).  Lrc potential
 * is based on integration of energy from cutoff to infinity, assuming no
 * pair correlations beyond the cutoff.
 */
public class P2SoftSphericalTruncatedShifted extends P2SoftSphericalTruncated {
    
    public P2SoftSphericalTruncatedShifted(Potential2SoftSpherical potential, double truncationRadius) {
        super(potential, truncationRadius);
    }
    
    /**
     * Returns the wrapped potential.
     */
    public Potential2SoftSpherical getWrappedPotential() {
        return potential;
    }
 
    /**
     * Returns the energy of the wrapped potential if the separation
     * is less than the cutoff value
     * @param r2 the squared distance between the atoms
     */
    public double u(double r2) {
        return (r2 < r2Cutoff) ? (potential.u(r2) - shift) : 0.0;
    }

    /**
     * Mutator method for the radial cutoff distance.
     */
    public void setTruncationRadius(double rCut) {
        super.setTruncationRadius(rCut);
        shift = potential.u(r2Cutoff);
    }

    /**
     * Returns null because the shift can't be corrected.
     */
    public Potential0Lrc makeLrcPotential(IAtomType[] types) {
        return null;
    }
    
    private static final long serialVersionUID = 1L;
    protected double shift;
}
