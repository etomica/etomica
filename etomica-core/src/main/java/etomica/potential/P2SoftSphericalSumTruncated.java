/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Space;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;


/**
 * Wraps a soft-spherical potential(s) to apply a truncation to it.  Energy and
 * its derivatives are set to zero at a specified cutoff.  (No accounting is
 * made of the infinite force existing at the cutoff point).  Lrc potential
 * is based on integration of energy from cutoff to infinity, assuming no
 * pair correlations beyond the cutoff.  This class allows P2SoftSphericalSum3
 * to do all the work with the wrapped potential(s).
 */
public class P2SoftSphericalSumTruncated extends P2SoftSphericalSum {

    protected double rCutoff, r2Cutoff;

    public P2SoftSphericalSumTruncated(Space _space, double truncationRadius, Potential2Soft... potential) {
        super(_space, potential);
        setTruncationRadius(truncationRadius);
    }

    /**
     * Returns the energy of the wrapped potential if the separation
     * is less than the cutoff value
     *
     * @param r2 the squared distance between the atoms
     */
    public double u(double r2) {
        if (r2 > r2Cutoff) return 0;
        return super.uWrapped(r2);
    }

    /**
     * Returns the derivative (r du/dr) of the wrapped potential if the separation
     * is less than the cutoff value
     * @param r2 the squared distance between the atoms
     */
    public double du(double r2) {
        if (r2 > r2Cutoff) return 0;
        return duWrapped(r2);
    }

    public void u012add(double r2, double[] u012) {
        if (r2 > r2Cutoff) return;
        super.uduWrapped(r2, u012);
    }

    /**
     * Returns the 2nd derivative (r^2 d^2u/dr^2) of the wrapped potential if the separation
     * is less than the cutoff value
     * @param r2 the squared distance between the atoms
     */
    public double d2u(double r2) {
        if (r2 > r2Cutoff) return 0;
        return super.d2uWrapped(r2);
    }

    public double integral(double rC) {
        throw new RuntimeException("nope");
    }

    /**
     * Accessor method for the radial cutoff distance.
     */
    public double getTruncationRadius() {
        return rCutoff;
    }

    /**
     * Mutator method for the radial cutoff distance.
     */
    public void setTruncationRadius(double rCut) {
        rCutoff = rCut;
        r2Cutoff = rCut*rCut;
    }

    /**
     * Returns the truncation radius.
     */
    public double getRange() {
        return rCutoff;
    }

    @Override
    public void u01TruncationCorrection(double[] uCorrection, double[] duCorrection) {
        double A = space.sphereArea(1.0);
        double D = space.D();
        double u = uWrapped(r2Cutoff);
        double integral = integralWrapped(rCutoff);
        uCorrection[0] = integral;
        duCorrection[0] = (-A * space.powerD(rCutoff) * u - D * integral);
    }

    /**
     * Returns the dimension (length) of the radial cutoff distance.
     */
    public Dimension getTruncationRadiusDimension() {
        return Length.DIMENSION;
    }
}
