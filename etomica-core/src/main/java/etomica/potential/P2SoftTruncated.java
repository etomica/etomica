/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.space.Space;
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
public class P2SoftTruncated implements Potential2Soft {

    protected final Space space;
    protected final Vector dr;
    protected final Potential2Soft wrappedPotential;
    protected double rCutoff, r2Cutoff;
    protected Boundary boundary;

    public P2SoftTruncated(Potential2Soft potential, double truncationRadius, Space _space) {
        this.space = _space;
        this.wrappedPotential = potential;
        setTruncationRadius(truncationRadius);
        dr = space.makeVector();
    }
    
    /**
     * Returns the wrapped potential.
     */
    public Potential2Soft getWrappedPotential() {
        return wrappedPotential;
    }

    /**
     * Returns the energy of the wrapped potential if the separation
     * is less than the cutoff value
     */
    public double energy(IAtomList atoms) {
        dr.Ev1Mv2(atoms.get(1).getPosition(), atoms.get(0).getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        return (r2 < r2Cutoff) ? wrappedPotential.energy(atoms) : 0;
    }

    @Override
    public double u(Vector dr12, IAtom atom1, IAtom atom2) {
        return dr12.squared() < r2Cutoff ? wrappedPotential.u(dr12, atom1, atom2) : 0;
    }

    @Override
    public double udu(Vector dr12, IAtom atom1, IAtom atom2, Vector f1, Vector f2) {
        if (dr12.squared() > r2Cutoff) return 0;
        return wrappedPotential.udu(dr12, atom1, atom2, f1, f2);
    }

    @Override
    public double uduTorque(Vector dr12, IAtom atom1, IAtom atom2, Vector f1, Vector f2, Vector t1, Vector t2) {
        if (dr12.squared() > r2Cutoff) return 0;
        return wrappedPotential.uduTorque(dr12, atom1, atom2, f1, f2, t1, t2);
    }

    /**
     * Returns the value of uInt for the wrapped potential.
     */
    public double integral(double rC) {
        return wrappedPotential.integral(rC);
    }

    public double u(double r2) {
        return wrappedPotential.u(r2);
    }

    public double du(double r2) {
        return wrappedPotential.du(r2);
    }

    public double d2u(double r2) {
        return wrappedPotential.d2u(r2);
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

    /**
     * Returns the dimension (length) of the radial cutoff distance.
     */
    public Dimension getTruncationRadiusDimension() {
        return Length.DIMENSION;
    }

    @Override
    public void u01TruncationCorrection(double[] uCorrection, double[] duCorrection) {
        double A = space.sphereArea(1.0);
        double D = space.D();
        double u = wrappedPotential.u(r2Cutoff);
        double integral = wrappedPotential.integral(rCutoff);
        uCorrection[0] = integral;
        duCorrection[0] = (-A * space.powerD(rCutoff) * u - D * integral);
    }
}
