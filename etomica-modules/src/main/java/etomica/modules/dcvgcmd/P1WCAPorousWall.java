/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.dcvgcmd;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.potential.IPotentialField;
import etomica.space.Boundary;
import etomica.space.Vector;

/**
 * This acts as a 1-body WCA potential wall perpendicular to the z direction
 * with a hole.  The size and placement of the hole can be set.  The interior
 * of the hole is discontinuous; an atom in the hole that attempts to move out
 * the hole and into the wall will experience a discontinuity.
 */
public class P1WCAPorousWall implements IPotentialField {

    private double sigma, sigma2;
    private double epsilon;
    private double cutoff, cutoff2;
    private double poreRadius, poreRadius2;
    private double z;
    private final Boundary boundary;

    public P1WCAPorousWall(Box box, double sigma, double epsilon, double poreRadius, double z) {
        boundary = box.getBoundary();
        setSigma(sigma);
        setEpsilon(epsilon);
        setPoreRadius(poreRadius);
        setZ(z);
    }

    public double u(IAtom atom) {
        Vector r = atom.getPosition();
        double rz = r.getX(2);
        double dz = rz - Math.signum(rz) * z;
        double L = boundary.getBoxSize().getX(2);
        double Ldz = rz - Math.signum(rz) * L / 2;
        double dz2;
        if (Math.abs(Ldz) < Math.abs(dz)) {
            dz2 = Ldz * Ldz;
        } else {
            if (inPore(r)) return 0;
            dz2 = dz * dz;
        }
        if (dz2 > cutoff2) return 0;
        return energy(dz2);
    }

    /**
     * Computes the force (and adds it to f) for IAtom atom and returns the
     * energy due to the field.
     */
    public double udu(IAtom atom, Vector f) {
        Vector r = atom.getPosition();
        double rz = r.getX(2);
        double dz = rz - z * Math.signum(rz);
        double L = boundary.getBoxSize().getX(2);
        double Ldz = rz - L / 2 * Math.signum(rz);
        double dz2;
        if (Math.abs(Ldz) < Math.abs(dz)) {
            dz = Ldz;
        } else {
            if (inPore(r)) return 0;
        }
        dz2 = dz * dz;

        if (dz2 > cutoff2) return 0.0;

        double s2 = sigma2 / dz2;
        double s6 = s2 * s2 * s2;
        double u = 4 * epsilon * s6 * (s6 - 1.0) + epsilon;

        double fz = (48 * epsilon * s6 * (s6 - 0.5)) / dz;

        f.setX(2, f.getX(2) + fz);

        return u;
    }

    private double energy(double r2) {
        r2 = sigma2 / r2;
        double r6 = r2 * r2 * r2;
        return 4 * epsilon * r6 * (r6 - 1.0) + epsilon;
    }

    private boolean inPore(Vector r) {
        double dx = r.getX(0);
        double dy = r.getX(1);
        double r2 = dx * dx + dy * dy;
        return r2 < poreRadius2;
    }

    /**
     * Returns the radius.
     *
     * @return double
     */
    public double getSigma() {
        return sigma;
    }

    /**
     * Sets the radius.
     *
     * @param radius The radius to set
     */
    public void setSigma(double radius) {
        this.sigma = radius;
        sigma2 = radius * radius;
        cutoff = radius * Math.pow(2, 1. / 6.);
        cutoff2 = cutoff * cutoff;
    }

    /**
     * @return Returns the epsilon.
     */
    public double getEpsilon() {
        return epsilon;
    }

    /**
     * @param epsilon The epsilon to set.
     */
    public void setEpsilon(double epsilon) {
        this.epsilon = epsilon;
    }


    /**
     * @return Returns the poreRadius.
     */
    public double getPoreRadius() {
        return poreRadius;
    }

    /**
     * @param poreRadius The poreRadius to set.
     */
    protected void setPoreRadius(double poreRadius) {
        this.poreRadius = poreRadius;
        poreRadius2 = poreRadius * poreRadius;
    }


    /**
     * @return Returns the z.
     */
    public double getZ() {
        return z;
    }

    /**
     * @param z The z to set.
     */
    protected void setZ(double z) {
        this.z = z;
    }
}
