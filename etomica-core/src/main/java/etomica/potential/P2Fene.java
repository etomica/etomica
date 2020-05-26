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
 * Finite elastic nonlinear extensible (FENE) spring potential.
 * Intramolecular potential for modeling polymer chains.
 *
 * @author David Kofke
 */
public class P2Fene extends Potential2SoftSpherical {

    private static final long serialVersionUID = 1L;
    private double r0, r02, h, prefactor;
    
    public P2Fene(Space _space) {
        this(_space, 1.50, 30.0);
    }
    
    public P2Fene(Space _space, double r0, double amplitude) {
        super(_space);
        setMaximumSeparation(r0);
        setAmplitude(amplitude);
    }

    public double u(double r2) {
        return (r2 < r02) ? prefactor * Math.log(1 - r2/r02) : Double.POSITIVE_INFINITY;
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        if (r2 > 0.99 * r02) return -h * r02 * 10;
        return h*r2*r02/(r02 - r2);
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
        double d = (r02 - r2);
        return h*r2*r02*(r02 + r2)/(d*d);
    }
            
    /**
     *  Integral used for corrections to potential truncation.
     */
    public double uInt(double rC) {
        return 0.0;
    }
    
    public double getAmplitude() {return h;}
    public void setAmplitude(double H) {
        this.h = H;
        prefactor = -0.5*h*r02;
    }
    
    /**
     * Not implemented correctly.  
     * Should be energy/length^3.
     */
    public Dimension getAmplitudeDimension() {
        return new CompoundDimension(new Dimension[]{Energy.DIMENSION,Length.DIMENSION},new double[]{1,-3});
    }

    /**
     * Accessor method for harmonic energy parameter
     */
    public double getMaximumSeparation() {return r0;}
    /**
     * Accessor method for harmonic energy parameter
     */
    public void setMaximumSeparation(double r0) {
        this.r0 = r0;
        r02 = r0*r0;
        prefactor = -0.5*h*r02;
    }
    
    public Dimension getMaximumSeparationDimension() {
        return Length.DIMENSION;
    }
    
}
