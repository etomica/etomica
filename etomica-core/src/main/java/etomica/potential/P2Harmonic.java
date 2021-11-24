/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Space;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;

/**
 * Harmonic Well interatomic potential.
 * Spherically symmetric potential of the form u(r) = 0.5*springConstant*(r-r0)^2 
 * where springConstant describes the strength of the pair interaction.
 *
 * @author Jhumpa Adhikari
 * @author David Kofke
 */
public class P2Harmonic extends Potential2SoftSpherical {

    private double w = 100.0;// Spring constant gives a measure of the strength of harmonic interaction
	private final boolean r0Zero;
	private double r0;
    
    public P2Harmonic(Space space, double w) {
    	this(space, w, 0.0);
    }
    /**
     * 
     * @param space
     * @param w spring constant
     * @param r0  Separation at which potential is at its minimum.  Default is
     * zero.
     */
    public P2Harmonic(Space space, double w, double r0) {
        super(space);
        setSpringConstant(w);
        r0Zero = (r0 == 0.0);
        setR0(r0);
    }

    public void u012add(double r2, double[] u012) {
        if (r0Zero) {
            u012[0] = u012[1] = u012[2] = w*r2;
            u012[0] *= 0.5;
            return;
        }
        double r = Math.sqrt(r2);
        double dx = r - r0;
        u012[0] = 0.5*w*dx*dx;
        u012[1] = w*r*dx;
        u012[2] = w*r2;
    }

    public double u(double r2) {
    	if(r0Zero) return 0.5*w*r2;
    	double dx = Math.sqrt(r2) - r0;
    	return 0.5*w*dx*dx;
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
    	if(r0Zero) return w*r2;
    	double r = Math.sqrt(r2);
    	return w*r*(r-r0);
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
        return w*r2;
    }
            
    /**
     * Integral used for corrections to potential truncation.
     */
    public double integral(double rC) {
        return 0.0;
    }

    /**
     * Accessor method for harmonic energy parameter
     */
    public double getSpringConstant() {return w;}
    /**
     * Accessor method for harmonic energy parameter
     */
    public void setSpringConstant(double factor) {
        w = factor;
    }

    /**
     * Not implemented correctly.  
     * Should be energy/length^2.
     */
    public Dimension getSpringConstantDimension() {
        return new CompoundDimension(new Dimension[]{Energy.DIMENSION,Length.DIMENSION},new double[]{1,-2});
    }
    
	/**
	 * Separation at which potential is at its minimum.
	 * @return double
	 */
	public double getR0() {
		return r0;
	}

	/**
	 * Sets the the separation at which potential is at its minimum.
	 * @param r0 The r0 to set
	 */
	public void setR0(double r0) {
		this.r0 = r0;
	}

    public Dimension getR0Dimension() {
        return Length.DIMENSION;
    }
    
}//end of P2Harmonic
