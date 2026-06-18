/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.exception.MethodNotImplementedException;
import etomica.space.Space;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;

/**
 * The Buckingham Exponential-6 atom-atom dispersion potential. Given formula:
 *
 * U(r) = epsilon*alpha/(alpha-6)*[(6/alpha)exp(alpha*[1-r/rmax]-(rmax/r)^6]
 * where epsilon describes the strength of the pair interaction, 
 * alpha is the repulsive steepness of the potential
 * rm is the distance which the potential is minimum 
 * and rmax is the point which the potential is maximum
 * @author Hye Min
 */

public class P2Exp6Buckingham implements IPotential2 {

    public static IPotential2 makeTruncated(double epsilon, double alpha, double rm, double rmax, TruncationFactory tf) {
        return tf.make(new P2Exp6Buckingham(epsilon, alpha, rm, rmax));
    }

    public P2Exp6Buckingham(double epsilon, double alpha, double rm, double rmax) {
        setEpsilon(epsilon);
        setAlpha(alpha);
        setRm(rm);
        setRmax(rmax);
    }

    public double getRmax() {
		return rmax;
	}

	public void setRmax(double rmax) {
		this.rmax = rmax;
	}

    /**
     * The energy u
     */
    public double u(double r2) {
        double r = Math.sqrt(r2);
        if (r < rmax) return Double.POSITIVE_INFINITY;
        double s = rm/r;
        double s2 = s * s;
        double s6 = s2*s2*s2;

        return  epsilon*alpha/(alpha-6)*(6/alpha* Math.exp(alpha*(1 - 1/s)) - s6);
        
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
    	double r = Math.sqrt(r2);
        if (r < rmax) return Double.POSITIVE_INFINITY;
        double s = rm/r;
        double s2 = s * s;
        double s7 = s2*s2*s2*s;
        return -6 * epsilon * alpha / s / (alpha - 6) * (Math.exp(alpha * (1 - 1 / s)) - s7);

    }

    /**
     * The second derivative of the pair energy, times the square of the
     * separation: r^2 d^2u/dr^2.
     */
    public double d2u(double r2) {
    	double r = Math.sqrt(r2);
        if (r < rmax) return Double.POSITIVE_INFINITY;
        double s = rm/r;
        double s2 = s * s;
        double s8 = s2*s2*s2*s2;
        return 6 * epsilon * alpha / s2 / (alpha - 6) * (alpha * Math.exp(alpha * (1 - 1 / s)) - 7 * s8);
    }

    /**
     * Integral used for corrections to potential truncation.
     */
    public double integral(Space space, double rC) { // need long range correction!!!!
        throw new MethodNotImplementedException("Integral for long-range correction for Exp-6 not yet implemented");
    }

    public double getEpsilon() {
        return epsilon;
    }

    public final void setEpsilon(double eps) {
        epsilon = eps;
    }

    public double getAlpha() {
        return alpha;
    }

    public final void setAlpha(double alp) {
        alpha = alp;
    }

    public double getRm() {
        return rm;
    }

    public final void setRm(double rm) {
    	this.rm = rm;
    }

	public Dimension getADimension() {
        return Energy.DIMENSION;
    }

    public Dimension getBDimension() {
        return Length.DIMENSION;
    }
  
    public Dimension getCDimension() {
        return new CompoundDimension(new Dimension[] {Energy.DIMENSION, Length.DIMENSION}, new double[] {1.0, 6.0});
    }

    private double rm;
    private double epsilon;
    private double alpha;
    private double rmax;
}
