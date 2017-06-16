/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomType;
import etomica.space.Space;

/**
 * Wraps a soft-spherical potential to apply a truncation to it.  Energy and
 * its derivatives are set to zero at a specified cutoff.  (No accounting is
 * made of the infinite force existing at the cutoff point).  Lrc potential
 * is based on integration of energy from cutoff to infinity, assuming no
 * pair correlations beyond the cutoff.
 */
public class P2SoftSphericalTruncatedShifted extends P2SoftSphericalTruncated {

    private static final long serialVersionUID = 1L;
    protected double shift;
 
    public P2SoftSphericalTruncatedShifted(Space _space, Potential2SoftSpherical potential,
                                           double truncationRadius) {
        super(_space, potential, truncationRadius);
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
    public Potential0Lrc makeLrcPotential(AtomType[] types) {
        return null;
    }
}
