/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Space;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;

/**
 * Soft-sphere interatomic potential.
 * Characterized by the pairwise-additive, spherically symmetric,
 * inverse-power potential of the form:
 * <p>
 * u(r) = epsilon*[sigma/ r]^n
 * <p>
 * where epsilon:  describes the strength of the pair interaction,
 * sigma  :  is the atom size parameter
 * n    :  is the degree of softness, s=1/n,
 * eg: hard sphere when s=0 or n->infinity
 * <p>
 * sigma/r  : sig_r
 * <p>
 * This class supports floating-point exponents, but is significantly slower
 * than P2SoftSphere (which requires an integer exponent)
 */
public class P2SoftSphereFloat implements IPotential2 {

    public static IPotential2 makeTruncated(double sigma, double epsilon, int n, TruncationFactory tf) {
        return tf.make(new P2SoftSphereFloat(sigma, epsilon, n));
    }

    public P2SoftSphereFloat(double sigma, double epsilon, double n) {
        setSigma(sigma);
        setEpsilon(epsilon);
        this.n = n;
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
        double s2 = sigma2 / r2;
        return epsilon * Math.pow(s2, 0.5 * n);
    }

    public void u012add(double r2, double[] u012) {
        double s2 = sigma2 / r2;
        double u = epsilon * Math.pow(s2, 0.5 * n);
        u012[0] += u;
        u012[1] += -n * u;
        u012[2] += (n + 1) * n * u;
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        return -n * u(r2);
    }

    /**
     * The second derivative of the pair energy, times the square of the
     * separation:  r^2 d^2u/dr^2.
     */
    public double d2u(double r2) {
        return n * (n + 1) * u(r2);
    }

    /**
     * Integral used for corrections to potential truncation.
     */
    public double integral(Space space, double rC) {
        double A = space.sphereArea(1.0);  //multiplier for differential surface element
        int D = space.D();                 //spatial dimension
        double rCD = space.powerD(rC);

        double sig_rCn = sigma / rC;
        return epsilon * A * rCD * Math.pow(sig_rCn, n) / (n - D);
    }

    /**
     * Accessor method for soft-sphere size parameter.
     */
    public double getSigma() {
        return sigma;
    }

    /**
     * Mutator method for soft sphere size parameter.
     * Does not adjust potential cutoff for change in sigma.
     */
    public final void setSigma(double sig) {
        sigma = sig;
        sigma2 = sig * sig;
    }

    public Dimension getSigmaDimension() {
        return Length.DIMENSION;
    }

    /**
     * Accessor method for soft-sphere energy parameter
     */
    public double getEpsilon() {
        return epsilon;
    }

    /**
     * Mutator method for soft-sphere energy parameter
     */
    public final void setEpsilon(double eps) {
        epsilon = eps;
    }

    public Dimension getEpsilonDimension() {
        return Energy.DIMENSION;
    }

    /**
     * Accessor method for soft-sphere softness parameter
     */
    public double getExponent() {
        return n;
    }

    private double sigma, sigma2;
    private double epsilon;
    private final double n;
}
