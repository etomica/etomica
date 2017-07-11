/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Vector;
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

public class P2Exp6Buckingham extends Potential2SoftSpherical {

    public P2Exp6Buckingham(Space _space) {
              this(_space, 1.0, 1.0, 1.0, 1.0);
    }

    public P2Exp6Buckingham(Space _space, double epsilon, double alpha, double rm, double rmax) {
        super(_space);
        dr01 = space.makeVector();
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
        double s = rm/r;
        double s2 = rmSquared/r2;
        double s6 = s2*s2*s2;

        if (r <
        	rmax) {
            return Double.POSITIVE_INFINITY;
        }
        //System.out.println(r+ " "+epsilon*alpha/(alpha-6)*(6/alpha* Math.exp(alpha*(1 - 1/s)) - s6));
        //System.out.println(s);
        return  epsilon*alpha/(alpha-6)*(6/alpha* Math.exp(alpha*(1 - 1/s)) - s6);
        
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
    	double r = Math.sqrt(r2);
        double s = rm/r;
        double s2 = rmSquared/r2;
        double s7 = s2*s2*s2*s;
        return -epsilon6*alpha /s/(alpha - 6)* (Math.exp(alpha*(1 - 1/s)) - s7);

    }

    /**
     * The second derivative of the pair energy, times the square of the
     * separation: r^2 d^2u/dr^2.
     */
    public double d2u(double r2) {
    	double r = Math.sqrt(r2);
        double s = rm/r;
        double s2 = rmSquared/r2;
        double s8 = s2*s2*s2*s2;
        return epsilon6*alpha/s2/(alpha - 6) * (alpha*Math.exp(alpha*(1 - 1/s)) - 7*s8);
    }

    /**
     * Integral used for corrections to potential truncation.
     */
    public double uInt(double rC) { // need long range correction!!!!
        throw new MethodNotImplementedException("Integral for long-range correction for Exp-6 not yet implemented");
    }

    public double getEpsilon() {
        return epsilon;
    }

    public final void setEpsilon(double eps) {
        epsilon = eps;
        epsilon6 = eps*6.0;
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
        rmSquared = rm*rm;
    }
  
    public double getRmSquared() {
		return rmSquared;
	}

	public void setRmSquared(double rmSquared) {
		this.rmSquared = rmSquared;
	}

	public double getEpsilon6() {
		return epsilon6;
	}

	public void setEpsilon6(double epsilon6) {
		this.epsilon6 = epsilon6;
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

    protected final Vector dr01;

    private static final long serialVersionUID = 1L;
    private double rm, rmSquared;
    private double epsilon;
    private double epsilon6;
    private double alpha;
    private double rmax;


}
