/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Space;


/**
 * Wraps a soft-spherical potential to apply a truncation to it.  Energy and
 * its derivatives are set to zero at a specified cutoff.  (No accounting is
 * made of the infinite force existing at the cutoff point).  A shift is
 * applied to both the energy and force  so that the they both go to 0
 * continuously at the cutoff.
 */
public class P2SoftSphericalSumTruncatedForceShifted extends P2SoftSphericalSumTruncated {

    protected double uShift, ufShift;

    public P2SoftSphericalSumTruncatedForceShifted(Space _space, double truncationRadius, Potential2SoftSpherical... potential) {
        super(_space, truncationRadius, potential);
    }

    public double u(double r2) {
        if (r2 > r2Cutoff) return 0;
        return uShift + Math.sqrt(r2) * ufShift + uWrapped(r2);
    }

    public double du(double r2) {
        if (r2 > r2Cutoff) return 0;
        return Math.sqrt(r2) * ufShift + duWrapped(r2);
    }

    /**
     * Mutator method for the radial cutoff distance.
     */
    public void setTruncationRadius(double rCut) {
        super.setTruncationRadius(rCut);
        uShift = -uWrapped(r2Cutoff);
        ufShift = -duWrapped(r2Cutoff) / rCutoff;
        uShift -= rCutoff * ufShift;
    }

    @Override
    public void u01TruncationCorrection(double[] uCorrection, double[] duCorrection) {
        double A = space.sphereArea(1.0);
        double D = space.D();
        double integral = integralWrapped(rCutoff);
        double u = uWrapped(r2Cutoff);
        double fsDU = -A * ufShift * space.powerD(rCutoff) * rCutoff / (D + 1);
        uCorrection[0] = integral - space.sphereVolume(rCutoff) * uShift + fsDU;
        duCorrection[0] = -A * space.powerD(rCutoff) * u - D * integral + fsDU;
    }
}
