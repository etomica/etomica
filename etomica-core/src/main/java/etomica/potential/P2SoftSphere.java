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
 *  inverse-power potential of the form:
 *  
 * 			     u(r) = epsilon*[sigma/ r]^n 
 *
 * where epsilon:  describes the strength of the pair interaction, 
 *       sigma  :  is the atom size parameter
 *         n    :  is the degree of softness, s=1/n, 
 *               eg: hard sphere when s=0 or n->infinity
 *
 *     sigma/r  : sig_r
 *
 *
 * @author Tai Boon Tan
 */
public class P2SoftSphere extends Potential2SoftSpherical {

    public static Potential2Soft makeTruncated(Space space, double sigma, double epsilon, int n, TruncationFactory tf) {
        return tf.make(new P2SoftSphereFloat(space, sigma, epsilon, n));
    }

    public P2SoftSphere(Space space) {
        this(space, 1.0, 1.0, 12);

    }

    public P2SoftSphere(Space space, double sigma, double epsilon, int n) {

        super(space);
        setSigma(sigma);
        setEpsilon(epsilon);
        this.n = n;
        evenN = n%2 == 0;
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
        double s2 = sigma2 / r2;
        double sig_rn = 1;

        if (n >= 1e6) {
            return r2 < sigma2 ? 0.0 : Double.POSITIVE_INFINITY;
        }

        if (n > 100) {
            sig_rn = Math.pow(s2, 0.5 * n);

        } else {
            double s6 = s2 * s2 * s2;
            switch (n) {
                case 4:
                    sig_rn = s2 * s2;
                    break;
                case 5:
                    sig_rn = s2 * s2 * Math.sqrt(s2);
                    break;
                case 6:
                    sig_rn = s6;
                    break;
                case 8:
                    sig_rn = s6 * s2;
                    break;
                case 9:
                    sig_rn = s6 * s2 * Math.sqrt(s2);
                    break;
                case 10:
                    sig_rn = s6 * s2 * s2;
                    break;
                case 12:
                    sig_rn = s6 * s6;
                    break;
                case 15:
                    sig_rn = s6 * s6 * s2 * Math.sqrt(s2);
                    break;
                case 16:
                    sig_rn = s6 * s6 * s2 * s2;
                    break;
                case 20:
                    sig_rn = s6 * s6 * s6 * s2;
                    break;
                case 24:
                    sig_rn = s6 * s6 * s6 * s6;
                    break;
                default:
                    if (n > 36) {
                        double s24 = s6 * s6 * s6 * s6;
                        for (int i = 0; i < n / 24; i++) {
                            sig_rn *= s24;
                        }
                    }
                    for (int i = 0; i < (n % 24) / 6; i++) {
                        sig_rn *= s6;
                    }
                    for (int i = 0; i < (n % 6) / 2; i++) {
                        sig_rn *= s2;
                    }
                    if (!evenN) {
                        sig_rn *= Math.sqrt(s2);
                    }
                    break;
            }
    	}
    	
    	return epsilon*sig_rn;
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        return -n*u(r2);
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
        return n*(n+1)*u(r2);
    }
            
    /**
     * Integral used for corrections to potential truncation.
     */
    public double integral(double rC) {
        double A = space.sphereArea(1.0);  //multiplier for differential surface element
        int D = space.D();                 //spatial dimension
        double rC3 = rC * rC * rC;

        return A * rC3 * u(rC * rC) / (n - D);  //complete LRC is obtained by multiplying by N1*N2/V
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
    public void setSigma(double sig) {
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
    public void setEpsilon(double eps) {
        epsilon = eps;
    }

    public Dimension getEpsilonDimension() {
        return Energy.DIMENSION;
    }

    /**
     * Accessor method for soft-sphere softness parameter
     */
    public int getExponent() {
        return n;
    }

    protected double sigma, sigma2;
    protected double epsilon;
    protected final int n;
    protected final boolean evenN;
}
