/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.space.Space;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;


/**
 * Wraps an array of soft-spherical potentials to apply a truncation.  Energy
 * and its derivatives are set to zero at a specified cutoff.  (No accounting
 * is made of the infinite force existing at the cutoff point).
 */
public class P2SoftSphericalTruncatedSum extends Potential2SoftSpherical
        implements PotentialTruncated {

    protected final Potential2SoftSpherical[] potentials;
    protected double rCutoff, r2Cutoff;

    public P2SoftSphericalTruncatedSum(Space _space, Potential2SoftSpherical[] potentials, double truncationRadius) {
        super(_space);
        this.potentials = potentials;
        setTruncationRadius(truncationRadius);
    }

    /**
     * Returns the wrapped potential.
     */
    public Potential2SoftSpherical[] getWrappedPotentials() {
        return potentials;
    }

    public void setBox(Box box) {
        throw new UnsupportedOperationException();
    }

    /**
     * Returns the energy of the wrapped potential if the separation
     * is less than the cutoff value
     *
     * @param r2 the squared distance between the atoms
     */
    public double u(double r2) {
        if (r2 > r2Cutoff) return 0;
        double sum = 0;
        for (Potential2SoftSpherical potential : potentials) {
            sum += potential.u(r2);
        }
        return sum;
    }

    /**
     * Returns the derivative (r du/dr) of the wrapped potential if the separation
     * is less than the cutoff value
     *
     * @param r2 the squared distance between the atoms
     */
    public double du(double r2) {
        if (r2 > r2Cutoff) return 0;
        double sum = 0;
        for (Potential2SoftSpherical potential : potentials) {
            sum += potential.du(r2);
        }
        return sum;
    }

    public void udu(double r2, double[] u, double[] du) {
        if (r2 > r2Cutoff) return;
        for (Potential2SoftSpherical potential : potentials) {
            potential.udu(r2, u, du);
        }
    }

    /**
     * Returns the 2nd derivative (r^2 d^2u/dr^2) of the wrapped potential if the separation
     * is less than the cutoff value
     *
     * @param r2 the squared distance between the atoms
     */
    public double d2u(double r2) {
        if (r2 > r2Cutoff) return 0;
        double sum = 0;
        for (Potential2SoftSpherical potential : potentials) {
            sum += potential.du(r2);
        }
        return sum;
    }

    /**
     * Returns the value of uInt for the wrapped potential.
     */
    public double uInt(double rC) {
        double sum = 0;
        for (Potential2SoftSpherical potential : potentials) {
            sum += potential.uInt(rC);
        }
        return sum;

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
        r2Cutoff = rCut * rCut;
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
        double uSum = 0, integral = 0;
        for (Potential2SoftSpherical potential : potentials) {
            uSum += potential.u(r2Cutoff);
            integral += potential.integral(r2Cutoff);
        }
        uCorrection[0] = integral;
        duCorrection[0] = (-A * space.powerD(rCutoff) * (uSum) - D * integral);
    }

    /**
     * Returns the dimension (length) of the radial cutoff distance.
     */
    public Dimension getTruncationRadiusDimension() {
        return Length.DIMENSION;
    }

    /**
     * Returns the zero-body potential that evaluates the contribution to the
     * energy and its derivatives from pairs that are separated by a distance
     * exceeding the truncation radius.
     */
    public Potential0Lrc makeLrcPotential(AtomType[] types) {
        throw new UnsupportedOperationException();
    }

    public boolean getMakeLrc() {
        throw new UnsupportedOperationException();
    }

    public void setMakeLrc(boolean newMakeLrc) {
        throw new UnsupportedOperationException();
    }
}
