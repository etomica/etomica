/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Vector;
import etomica.space.Space;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Length;

/**
 * Dreiding: Lennard-Jones non-bonding potential.
 * Spherically symmetric potential of the form u(r) = D0*[(r0/r)^12 - 2*(r0/r)^6]
 * where D0 describes the van der Waals well depth [unit Kelvin ], 
 * and r0 is the van der Waals bond length [unit Amstrom ].
 * 
 * r0 = sigma; D0 = epsilon
 *
 * @author Tai Tan
 */

public class P2LennardJonesDreiding extends Potential2SoftSpherical {
	
	public P2LennardJonesDreiding(Space space) {
        this(space, 1.0, 1.0);
    }
	
    public P2LennardJonesDreiding(Space space, double sigma, double epsilon) {
        super(space);
        dr01 = space.makeVector();
        setSigma(sigma);
        setEpsilon(epsilon);
    
    }

    /**
     * The energy u.
     * rho = sigma / r
     */
    public double u(double r2) {
        double s2 = sigmaSquared/r2;
        rho6 = s2*s2*s2;
        return epsilon*rho6*(rho6 - 2.0) ;
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        double s2 = sigmaSquared/r2;
        rho6 = s2*s2*s2;
        return -epsilon12*rho6*(rho6 - 1.0);
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
        double s2 = sigmaSquared/r2;
        rho6 = s2*s2*s2;
        return epsilon156*rho6*(rho6 - _84div156);
    }
    
    /**
     *  Integral used for corrections to potential truncation.
     */
    public double uInt(double rC) { //need long range correction!!!!
        return 0.0;
    }

    /**
     * Accessor method for Lennard-Jones size parameter.
     */
    public double getSigma() {return sigma;}
    /**
     * Mutator method for Lennard-Jones size parameter.
     * Does not adjust potential cutoff for change in sigma.
     */
    public final void setSigma(double s) {
        sigma = s;
        sigmaSquared = s*s;
    }
    
    /**
     * Accessor method for Lennard-Jones energy parameter
     */
    public double getEpsilon() {return epsilon;}
    /**
     * Mutator method for Lennard-Jones energy parameter
     */
    public final void setEpsilon(double eps) {
        epsilon = eps;
        epsilon12 = eps*12.0;
        epsilon156 = eps*156.0;
        }
    public Dimension getSigmaDimension() {return Length.DIMENSION;}
    public Dimension getEpsilonDimension() {return Energy.DIMENSION;}
   
    private double sigma, sigmaSquared;
    private double epsilon;
    private double epsilon12;
    private double epsilon156;
    private static final double _84div156 = 84./156.;
    private double rho6;
    protected final Vector dr01;
	
	private static final long serialVersionUID = 1L;
}
