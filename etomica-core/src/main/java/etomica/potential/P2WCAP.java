/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;

/**
 * LJ potential minus WCA.
 */
public class P2WCAP implements Potential2Soft {

    protected Boundary boundary;
    protected final Vector dr;

    /**
     * Constructs potential using default sigma and epsilon given by Default class.
     */
    public P2WCAP(Space space) {
        this(space, 1.0, 1.0);
    }

    public P2WCAP(Space space, double sigma, double epsilon) {
        dr = space.makeVector();
        setSigma(sigma);
        setEpsilon(epsilon);
    }

    /**
     * Energy of the pair as given by the u(double) method
     */
    public double energy(IAtomList atoms) {
        dr.Ev1Mv2(atoms.get(1).getPosition(), atoms.get(0).getPosition());
        boundary.nearestImage(dr);
        return u(dr.squared());
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
        if (r2 < rangeSquared) {
            return -epsilon;
        }

        double s2 = sigmaSquared / r2;
        double s6 = s2 * s2 * s2;
//        System.out.println("uP: "+(epsilon4*s6*(s6-1)));
        return epsilon4 * s6 * (s6 - 1.0);
    }

    /**
     * Accessor method for the size parameter.
     */
    public double getSigma() {
        return sigma;
    }

    /**
     * Mutator method for Lennard-Jones size parameter.
     * Does not adjust potential cutoff for change in sigma.
     */
    public final void setSigma(double s) {
        sigma = s;
        sigmaSquared = s * s;
        range = sigma * Math.pow(2, 1. / 6.);
        rangeSquared = range * range;

    }

    public Dimension getSigmaDimension() {
        return Length.DIMENSION;
    }

    /**
     * Accessor method for the energy parameter
     */
    public double getEpsilon() {
        return epsilon;
    }

    /**
     * Mutator method for the energy parameter
     */
    public final void setEpsilon(double eps) {
        epsilon = eps;
        epsilon4 = eps * 4.0;
    }

    public Dimension getEpsilonDimension() {
        return Energy.DIMENSION;
    }

    private double sigma, sigmaSquared, range, rangeSquared;
    private double epsilon, epsilon4;
}
