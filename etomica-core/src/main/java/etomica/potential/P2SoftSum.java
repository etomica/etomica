/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Space;
import etomica.space.Vector;


/**
 * Wraps up to 3 soft-spherical potential and sums them.  This class can
 * happily take a 1 or 2 potentials and use them at no additional cost
 * (compared to a class that is hardcoded to take the correct number).  The
 * overhead for this class is much smaller than the overhead if holding and
 * iterating over an array of potentials.
 */
public class P2SoftSum implements Potential2Soft {

    private final Potential2Soft potential1, potential2, potential3;

    public P2SoftSum(Potential2Soft... potential) {
        if (potential.length > 3) throw new RuntimeException("This class only handles up to 3 potentials");
        this.potential1 = potential[0];
        this.potential2 = potential.length > 1 ? potential[1] : null;
        this.potential3 = potential.length > 2 ? potential[2] : null;
    }

    /**
     * Returns the wrapped potential.
     */
    public Potential2Soft getWrappedPotential1() {
        return potential1;
    }

    public Potential2Soft getWrappedPotential2() {
        return potential2;
    }

    public Potential2Soft getWrappedPotential3() {
        return potential3;
    }

    public double u(double r2) {
        return uWrapped(r2);
    }

    protected double uWrapped(double r2) {
        double u = potential1.u(r2);
        if (potential2 == null) return u;
        u += potential2.u(r2);
        if (potential3 == null) return u;
        return u + potential3.u(r2);
    }

    public double du(double r2) {
        return duWrapped(r2);
    }

    protected double duWrapped(double r2) {
        double du = potential1.du(r2);
        if (potential2 == null) return du;
        du += potential2.du(r2);
        if (potential3 == null) return du;
        return du + potential3.du(r2);
    }

    public void u012add(double r2, double[] u012) {
        uduWrapped(r2, u012);
    }

    protected void uduWrapped(double r2, double[] u012) {
        potential1.u012add(r2, u012);
        if (potential2 == null) return;
        potential2.u012add(r2, u012);
        if (potential3 == null) return;
        potential3.u012add(r2, u012);
    }

    public double d2u(double r2) {
        return d2uWrapped(r2);
    }

    protected double uWrapped(Vector dr12, IAtom atom1, IAtom atom2) {
        double u = potential1.u(dr12, atom1, atom2);
        if (potential2 == null) return u;
        u += potential2.u(dr12, atom1, atom2);
        if (potential3 == null) return u;
        return u + potential3.u(dr12, atom1, atom2);
    }

    @Override
    public double u(Vector dr12, IAtom atom1, IAtom atom2) {
        return uWrapped(dr12, atom1, atom2);
    }

    protected double uduWrapped(Vector dr12, IAtom atom1, IAtom atom2, Vector f1, Vector f2) {
        double u = potential1.udu(dr12, atom1, atom2, f1, f2);
        if (potential2 == null) return u;
        u += potential2.udu(dr12, atom1, atom2, f1, f2);
        if (potential3 == null) return u;
        return u + potential3.udu(dr12, atom1, atom2, f1, f2);
    }

    @Override
    public double udu(Vector dr12, IAtom atom1, IAtom atom2, Vector f1, Vector f2) {
        return uduWrapped(dr12, atom1, atom2, f1, f2);
    }

    protected double uduTorqueWrapped(Vector dr12, IAtom atom1, IAtom atom2, Vector f1, Vector f2, Vector t1, Vector t2) {
        double u = potential1.uduTorque(dr12, atom1, atom2, f1, f2, t1, t2);
        if (potential2 == null) return u;
        u += potential2.uduTorque(dr12, atom1, atom2, f1, f2, t1, t2);
        if (potential3 == null) return u;
        return u + potential3.uduTorque(dr12, atom1, atom2, f1, f2, t1, t2);
    }

    @Override
    public double uduTorque(Vector dr12, IAtom atom1, IAtom atom2, Vector f1, Vector f2, Vector t1, Vector t2) {
        return uduTorqueWrapped(dr12, atom1, atom2, f1, f2, t1, t2);
    }

    protected Hessian d2uWrapped(Vector dr12, IAtom atom1, IAtom atom2) {
        Hessian h1 = potential1.d2u(dr12, atom1, atom2);
        if (potential2 == null) return h1;
        Hessian h2 = potential2.d2u(dr12, atom1, atom2);
        if (h2 != null) {
            if (h1 == null) {
                h1 = h2;
            }
            else {
                h1.PE(h2);
            }
        }
        if (potential3 == null) return h1;
        Hessian h3 = potential3.d2u(dr12, atom1, atom2);
        if (h3 != null) {
            if (h1 == null) {
                h1 = h3;
            } else {
                h1.PE(h3);
            }
        }
        return h1;
    }

    @Override
    public Hessian d2u(Vector dr12, IAtom atom1, IAtom atom2) {
        return d2uWrapped(dr12, atom1, atom2);
    }

    protected double d2uWrapped(double r2) {
        double d2u = potential1.d2u(r2);
        if (potential2 == null) return d2u;
        d2u += potential2.d2u(r2);
        if (potential3 == null) return d2u;
        return d2u + potential3.d2u(r2);
    }

    public double integral(Space space, double rC) {
        return integralWrapped(space, rC);
    }

    public double integralWrapped(Space space, double rC) {
        double d2u = potential1.integral(space, rC);
        if (potential2 == null) return d2u;
        d2u += potential2.integral(space, rC);
        if (potential3 == null) return d2u;
        return d2u + potential3.integral(space, rC);
    }

    /**
     * Returns the truncation radius.
     */
    public double getRange() {
        double r = potential1.getRange();
        if (potential2 == null) return r;
        r = Math.max(r, potential2.getRange());
        if (potential3 == null) return r;
        return Math.max(r, potential3.getRange());
    }

    @Override
    public double energy(IAtomList atoms) {
        return 0;
    }

}
