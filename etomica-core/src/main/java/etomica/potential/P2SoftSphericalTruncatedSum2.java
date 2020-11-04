/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomType;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;


/**
 * Wraps a soft-spherical potential to apply a truncation to it.  Energy and
 * its derivatives are set to zero at a specified cutoff.  (No accounting is
 * made of the infinite force existing at the cutoff point).  Lrc potential
 * is based on integration of energy from cutoff to infinity, assuming no
 * pair correlations beyond the cutoff.
 */
public class P2SoftSphericalTruncatedSum2 extends Potential2SoftSpherical
               implements PotentialTruncated {

    protected final Potential2SoftSpherical potential1;
    protected final Potential2SoftSpherical potential2;
    protected double rCutoff, r2Cutoff;

    public P2SoftSphericalTruncatedSum2(Space _space, Potential2SoftSpherical potential1, Potential2SoftSpherical potential2, double truncationRadius) {
        super(_space);
        this.potential1 = potential1;
        this.potential2 = potential2;
        setTruncationRadius(truncationRadius);
    }

    /**
     * Returns the wrapped potential.
     */
    public Potential2SoftSpherical getWrappedPotential1() {
        return potential1;
    }

    public Potential2SoftSpherical getWrappedPotential2() {
        return potential2;
    }

    public void setBox(Box box) {
        throw new UnsupportedOperationException();
    }

    /**
     * Returns the energy of the wrapped potential if the separation
     * is less than the cutoff value
     * @param r2 the squared distance between the atoms
     */
    public double u(double r2) {
        return (r2 < r2Cutoff) ? potential1.u(r2) + potential2.u(r2) : 0.0;
    }

    /**
     * Returns the derivative (r du/dr) of the wrapped potential if the separation
     * is less than the cutoff value
     * @param r2 the squared distance between the atoms
     */
    public double du(double r2) {
        return (r2 < r2Cutoff) ? potential1.du(r2) + potential2.du(r2) : 0.0;
    }

    public void udu(double r2, double[] u, double[] du) {
        if (r2 > r2Cutoff) return;
        potential1.udu(r2, u, du);
        potential2.udu(r2, u, du);
    }

    /**
     * Returns the 2nd derivative (r^2 d^2u/dr^2) of the wrapped potential if the separation
     * is less than the cutoff value
     * @param r2 the squared distance between the atoms
     */
    public double d2u(double r2) {
        return (r2 < r2Cutoff) ? potential1.d2u(r2) + potential2.d2u(r2) : 0.0;
    }

    /**
     * Returns the value of uInt for the wrapped potential.
     */
    public double uInt(double rC) {
        return potential1.uInt(rC) + potential2.uInt(rC);
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
        double integral = potential1.integral(rCutoff) + potential2.integral(rCutoff);
        uCorrection[0] = integral;
        duCorrection[0] = (-A * space.powerD(rCutoff) * (potential1.u(r2Cutoff) + potential2.u(r2Cutoff)) - D * integral);
    }

    /**
     * Returns the dimension (length) of the radial cutoff distance.
     */
    public Dimension getTruncationRadiusDimension() {return Length.DIMENSION;}
    
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
