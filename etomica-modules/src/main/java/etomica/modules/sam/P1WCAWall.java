/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.sam;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.potential.IPotentialField;
import etomica.potential.Potential1;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * 1-D potential that has a WCA form in the Z direction.
 */

public class P1WCAWall extends Potential1 implements IPotentialField {

    protected double sigma, sigma2;
    protected double epsilon;
    protected double cutoff, cutoff2;
    protected int wallDim;
    protected double wallPosition;

    public P1WCAWall(Space space, int wallDim, double sigma, double epsilon) {
        super(space);
        setSigma(sigma);
        setEpsilon(epsilon);
        setWallDim(wallDim);
    }

    @Override
    public double u(IAtom atom) {
        double rz = atom.getPosition().getX(wallDim) - wallPosition;
        return energy(rz*rz);
    }

    public double getRange() {
        return cutoff;
    }

    public double energy(IAtomList atom) {
        double rz = atom.get(0).getPosition().getX(wallDim) - wallPosition;
        return energy(rz*rz);
    }

    private double energy(double r2) {
        if (r2 > cutoff2) {
            return 0;
        }
        double s2 = sigma2 / r2;
        double s6 = s2 * s2 * s2;
        return 4 * epsilon * s6 * (s6 - 1.0) + epsilon;
    }

    public double udu(IAtom atom, Vector force) {
        double rz = atom.getPosition().getX(wallDim) - wallPosition;
        double r2 = rz*rz;
        if (r2 > cutoff2) {
            return 0;
        }
        double s2 = sigma2 / r2;
        double s6 = s2 * s2 * s2;
        double u = 4 * epsilon * s6 * (s6 - 1.0) + epsilon;
        double gradz = -48 * epsilon * s6 * (s6 - 0.5);
        double fx = force.getX(0);
        fx += rz > 0 ? -gradz : +gradz;
        force.setX(0, fx);
        return u;
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
     * @param radius
     *            The radius to set
     */
    public void setSigma(double radius) {
        this.sigma = radius;
        sigma2 = sigma*sigma;
        cutoff2 = sigma2 * Math.pow(2, 1. / 3.);
    }

    /**
     * @return Returns the epsilon.
     */
    public double getEpsilon() {
        return epsilon;
    }

    /**
     * @param epsilon
     *            The epsilon to set.
     */
    public void setEpsilon(double epsilon) {
        this.epsilon = epsilon;
    }

    public int getWallDim() {
        return wallDim;
    }

    public void setWallDim(int wallDim) {
        this.wallDim = wallDim;
    }

    public double getWallPosition() {
        return wallPosition;
    }

    public void setWallPosition(double newWallPosition) {
        wallPosition = newWallPosition;
    }
}
