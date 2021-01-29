/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.space.Space;


/**
 * Wraps up to 3 soft-spherical potential and sums them.  This class can
 * happily take a 1 or 2 potentials and use them at no additional cost
 * (compared to a class that is hardcoded to take the correct number).  The
 * overhead for this class is much smaller than the overhead if holding and
 * iterating over an array of potentials.
 */
public class P2SoftSphericalSum extends Potential2SoftSpherical
        implements PotentialTruncated {

    private final Potential2SoftSpherical potential1, potential2, potential3;

    public P2SoftSphericalSum(Space _space, Potential2SoftSpherical[] potential) {
        this(_space, potential[0], potential.length > 1 ? potential[1] : null, potential.length > 2 ? potential[2] : null);
    }

    public P2SoftSphericalSum(Space _space, Potential2SoftSpherical potential1) {
        this(_space, potential1, null, null);
    }

    public P2SoftSphericalSum(Space _space, Potential2SoftSpherical potential1, Potential2SoftSpherical potential2) {
        this(_space, potential1, potential2, null);
    }

    public P2SoftSphericalSum(Space _space, Potential2SoftSpherical potential1, Potential2SoftSpherical potential2, Potential2SoftSpherical potential3) {
        super(_space);
        this.potential1 = potential1;
        this.potential2 = potential2;
        this.potential3 = potential3;
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

    public Potential2SoftSpherical getWrappedPotential3() {
        return potential3;
    }

    public void setBox(Box box) {
        throw new UnsupportedOperationException();
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

    public void udu(double r2, double[] u, double[] du) {
        uduWrapped(r2, u, du);
    }

    protected void uduWrapped(double r2, double[] u, double[] du) {
        potential1.udu(r2, u, du);
        if (potential2 == null) return;
        potential2.udu(r2, u, du);
        if (potential3 == null) return;
        potential3.udu(r2, u, du);
    }

    public double d2u(double r2) {
        return d2uWrapped(r2);
    }

    protected double d2uWrapped(double r2) {
        double d2u = potential1.d2u(r2);
        if (potential2 == null) return d2u;
        d2u += potential2.d2u(r2);
        if (potential3 == null) return d2u;
        return d2u + potential3.d2u(r2);
    }

    public double uInt(double rC) {
        return uIntWrapped(rC);
    }

    public double uIntWrapped(double rC) {
        double d2u = potential1.uInt(rC);
        if (potential2 == null) return d2u;
        d2u += potential2.uInt(rC);
        if (potential3 == null) return d2u;
        return d2u + potential3.uInt(rC);
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

    /**
     * Returns the zero-body potential that evaluates the contribution to the
     * energy and its derivatives from pairs that are separated by a distance
     * exceeding the truncation radius.
     */
    public Potential0Lrc makeLrcPotential(AtomType[] types) {
        throw new UnsupportedOperationException();
    }
}
