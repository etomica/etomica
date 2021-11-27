/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

/**
 * Soft-spherical potential class that shifts both the potential energy and the
 * force such that both are 0 at the cutoff.  The potential is of the form
 * U = U_LJ + Ar + B
 * 
 * @author Andrew Schultz
 */
public class P2SoftSphericalTruncatedForceShifted extends
        P2SoftSphericalTruncatedShifted {

    public P2SoftSphericalTruncatedForceShifted(IPotential2 potential, double truncationRadius) {
        super(potential, truncationRadius);
    }

    /**
     * Mutator method for the radial cutoff distance.
     */
    public void setTruncationRadius(double rCut) {
        super.setTruncationRadius(rCut);
        fShift = potential.du(r2Cutoff)/rCut;
        shift = potential.u(r2Cutoff) - fShift*rCut;
    }
    
    public double u(double r2) {
        return (r2 < r2Cutoff) ? (potential.u(r2) - fShift*Math.sqrt(r2) - shift) : 0.0;
    }
    
    public double du(double r2) {
        return (r2 < r2Cutoff) ? (potential.du(r2) - fShift*Math.sqrt(r2)) : 0.0;
    }

    public void u012add(double r2, double[] u012) {
        if (r2 > r2Cutoff) return;
        potential.u012add(r2, u012);
        double r = Math.sqrt(r2);
        u012[0] -= fShift * r + shift;
        u012[1] -= fShift * r;
    }

    protected double fShift;
}
