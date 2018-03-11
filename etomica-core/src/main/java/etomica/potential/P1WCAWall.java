/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * 1-D potential that has a WCA form in the Z direction.
 */

public class P1WCAWall extends Potential1 implements PotentialSoft {

    private static final long serialVersionUID = 1L;
    protected final Vector[] gradient;
    protected double sigma;
    protected double epsilon;
    protected double cutoff;
    protected int wallDim;

    public P1WCAWall(Space space, int wallDim) {
        this(space, wallDim, 1.0, 1.0);
    }

    public P1WCAWall(Space space, int wallDim, double sigma, double epsilon) {
        super(space);
        setSigma(sigma);
        setEpsilon(epsilon);
        setWallDim(wallDim);
        gradient = new Vector[1];
        gradient[0] = space.makeVector();
    }

    public double getRange() {
        return cutoff;
    }

    public double energy(IAtomList atom) {
        Vector dimensions = boundary.getBoxSize();
        double rz = atom.get(0).getPosition().getX(wallDim);
        double dzHalf = 0.5 * dimensions.getX(wallDim);
        return energy(dzHalf + rz) + energy(dzHalf - rz);
    }

    private double energy(double r) {
        if (r > cutoff) {
            return 0;
        }
        double rr = sigma / r;
        double r2 = rr * rr;
        double r6 = r2 * r2 * r2;
        return 4 * epsilon * r6 * (r6 - 1.0) + epsilon;
    }

    private double gradient(double r) {
        if (r > cutoff) {
            return 0;
        }
        double rr = sigma / r;
        double r2 = rr * rr;
        double r6 = r2 * r2 * r2;
        return -48 * epsilon * r6 * (r6 - 0.5);
    }

    public Vector[] gradient(IAtomList atom) {
        Vector dimensions = boundary.getBoxSize();
        double rz = atom.get(0).getPosition().getX(wallDim);
        double dzHalf = 0.5 * dimensions.getX(wallDim);
        double gradz = gradient(rz + dzHalf) - gradient(dzHalf - rz);
        gradient[0].setX(wallDim, gradz);
        return gradient;
    }
    
    public Vector[] gradient(IAtomList atom, Tensor pressureTensor) {
        return gradient(atom);
    }
    
    public double virial(IAtomList atoms) {
        return 0.0;
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
        cutoff = radius * Math.pow(2, 1. / 6.);
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
}
