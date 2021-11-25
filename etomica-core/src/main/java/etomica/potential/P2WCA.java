/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Space;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;

/**
 * Weeks-Chandler-Andersen potential.  Obtained by truncating the Lennard-Jones
 * potential at the separation where it has its minimum, and shifting it upwards so the
 * minimum is at zero energy.  The resulting potential has a value and first
 * derivative that is continuous at the point of truncation.
 *
 * @author David Kofke (edited by Eric C. Cichowski and Todd Schmidt)
 */
public class P2WCA implements Potential2Soft {

    /**
     * Constructs potential using default sigma and epsilon given by Default class.
     */
    public P2WCA() {
        this(1.0, 1.0);
    }
    
    public P2WCA(double sigma, double epsilon) {
        setSigma(sigma);
        setEpsilon(epsilon);
    }

    /**
     * Returns the range of the potential, which is the point of truncation.  This
     * is equal to 2^(1/6) * sigma.
     */
    public double getRange() {
        return range;
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
        if(r2 < rangeSquared) {
            double s2 = sigmaSquared/r2;
            s6 = s2*s2*s2;
            return epsilon4*s6*(s6 - 1.0) + epsilon;
        }
        return 0.0;
     }


    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        double du = 0.0;
        if(r2 < rangeSquared) {
            double s2 = sigmaSquared/r2;
            s6 = s2*s2*s2;
            du = -epsilon48*s6*(s6 - 0.5);
        }
//        System.out.println("r2, du: "+r2 + "  "+du);
        return du;
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
        if(r2 < rangeSquared) {
            double s2 = sigmaSquared/r2;
            s6 = s2*s2*s2;
            return epsilon624*s6*(s6 - _168div624);
        }
        return 0.0;
      }
            
    /**
     * Returns zero, since there is no long-range correction.
     */
    public double integral(Space space, double rC) {
        return 0.0;
    }

    /**
     * Accessor method for the size parameter.
     */
    public double getSigma() {return sigma;}
    /**
     * Mutator method for Lennard-Jones size parameter.
     * Does not adjust potential cutoff for change in sigma.
     */
    public final void setSigma(double s) {
        sigma = s;
        sigmaSquared = s*s;
		range = sigma*Math.pow(2,1./6.);
		rangeSquared = range*range;

    }
    public Dimension getSigmaDimension() {return Length.DIMENSION;}
    
    /**
     * Accessor method for the energy parameter
     */
    public double getEpsilon() {return epsilon;}
    
    /**
     * Mutator method for the energy parameter
     */
    public final void setEpsilon(double eps) {
        epsilon = eps;
        epsilon4 = eps*4.0;
        epsilon48 = eps*48.0;
        epsilon624 = eps*624.0;
    }
    public Dimension getEpsilonDimension() {return Energy.DIMENSION;}
   
    private double sigma, sigmaSquared, range, rangeSquared;
    private double epsilon;
    private double epsilon4, epsilon48, epsilon624;
    private static final double _168div624 = 168./624.;
    private double s6;
}
