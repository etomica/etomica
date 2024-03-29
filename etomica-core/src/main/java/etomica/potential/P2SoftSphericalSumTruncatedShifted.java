/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Space;


/**
 * Wraps a soft-spherical potential to apply a truncation to it.  Energy and
 * its derivatives are set to zero at a specified cutoff.  (No accounting is
 * made of the infinite force existing at the cutoff point).  A shift is
 * applied so that the energy goes to 0 continuously at the cutoff.
 */
public class P2SoftSphericalSumTruncatedShifted extends P2SoftSphericalSumTruncated {

    protected double shift;

    public P2SoftSphericalSumTruncatedShifted(double truncationRadius, IPotential2... potential) {
        super(truncationRadius, potential);
    }

    public double u(double r2) {
        if (r2 > r2Cutoff) return 0;
        return shift + uWrapped(r2);
    }

    public void u012add(double r2, double[] u012) {
        if (r2 > r2Cutoff) return;
        super.uduWrapped(r2, u012);
        u012[0] += shift;
    }

    /**
     * Mutator method for the radial cutoff distance.
     */
    public void setTruncationRadius(double rCut) {
        super.setTruncationRadius(rCut);
        shift = -uWrapped(r2Cutoff);
    }

    @Override
    public void u01TruncationCorrection(Space space, double[] uCorrection, double[] duCorrection) {
        double A = space.sphereArea(1.0);
        double D = space.D();
        double u = uWrapped(r2Cutoff);
        double integral = integralWrapped(space, rCutoff);
        uCorrection[0] = integral - space.sphereVolume(rCutoff) * shift;
        duCorrection[0] = -A * space.powerD(rCutoff) * u - D * integral;
    }
}
