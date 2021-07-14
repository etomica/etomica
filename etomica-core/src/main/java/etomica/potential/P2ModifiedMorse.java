/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;
import etomica.space.Space;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;
import etomica.util.Constants;

/*
 * Modified Morse potential
 * 
 *  The functional form of potential is :
 *  
 *  	u(i,j) = -epsilon(i,j) * ( 1 - { 1 - exp( -A(i,j) ( r(i,j) - re(i,j) ) }^2 ) + z(i)*z(j)*(electron_charge)^2/(4*pi*free_space_permittivity*r(i,j))	
 *   
 * where epsilon: describes the strength of the pair interaction - the well depth (energy) 
 *         re   : is the equilibrium pair separation - the position of the well bottom (distance)
 *         A    : is a parameter controlling the width of the potential well (inverse distance)          
 *  
 *
 * adapted from P2Morse.java by K.R. Schadel (2008)
 */


public final class P2ModifiedMorse extends Potential2SoftSpherical {

    public static Potential2Soft makeTruncated(Space space, double epsilon, double re, double a, double z1, double z2, double coulombicCutoff, TruncationFactory tf) {
        return tf.make(new P2ModifiedMorse(space, epsilon, re, a, z1, z2, coulombicCutoff));
    }

    public P2ModifiedMorse(Space space) {
        this(space, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    }

    public P2ModifiedMorse(Space space, double epsilon, double re, double a, double z1, double z2, double coulombicCutoff) {
        super(space);
        setEpsilon(epsilon);
        setRe(re);
        setA(a);
        setZ1(z1);
        setZ2(z2);
        this.coulombicCutoff = coulombicCutoff;
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
    	double r = Math.sqrt(r2);
    	
    	/* The O-aC and O-aH Coulombic interactions are highly attractive for small r.
    	 * As r goes to zero, the Coulombic interaction energy goes to negative infinity.
    	 * This causes some overlapped configurations to have aphysical attractive interaction energies.
    	 */
    	
    	
    	double expTerm = Math.exp(a*(re-r));
    	double morseTerm = epsilon*(expTerm-1)*(expTerm-1)-epsilon;
    	
    	
    	
    	// Note: e is the unit of charge in simulation units (electron charge rather than Coulombs)
    	double chargeTerm = z1*z2/(4*Math.PI*Constants.EPSILON_0*r);
    	
    	if (r <= coulombicCutoff) {
    		
    		return Double.POSITIVE_INFINITY;
    		
    		/*if (a == 5.42872) {
        		if (r > 1000) {
        			System.out.println("aC-aC distance = " + r + "Angstroms");
        		}
        	}*/
    	/*} else 	if (r < 100) {
    			return morseTerm + chargeTerm;
    	} else {
    		return morseTerm;
    		
        } 
    	*/
    		
	    }
    	return morseTerm + chargeTerm;
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
    	double alpha = a*re;
    	double r = Math.sqrt(r2);
    	double expTerm = Math.exp(-alpha*(-1+(r/re)));
    	
    	System.out.println("Ahh!!!!  This has not been modified yet!!!");
    	return -2*(r/re)*epsilon*alpha*(expTerm-1)*expTerm;
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
    	double alpha = a*re;
    	double r = Math.sqrt(r2);
    	double expTerm = Math.exp(-alpha*(-1+(r/re)));
    	
    	System.out.println("Ahh!!!!  This has not been modified yet!!!");
        return 2*(r2/(re*re))*epsilon*alpha*alpha*expTerm*(2*expTerm-1);
    }
            
    /**
     * Integral used for corrections to potential truncation.
     */
    public double integral(double rC) {
        double alpha = a * re;
        double A = space.sphereArea(1.0);  //multiplier for differential surface element
        double rC2 = rC * rC;
        double re2 = re * re;
        double alpha2 = alpha * alpha;
        double alpha3 = alpha * alpha * alpha;
        double expTerm = Math.exp(-alpha * (-1 + (rC / re)));

        System.out.println("Ahh!!!!  This has not been modified yet!!!");
        return (-A / (4 * alpha3))
        		*(re*expTerm*epsilon*
        				(2*expTerm*alpha2*rC2 + 2*expTerm*alpha*re*rC	
        						+expTerm*re2 -8*alpha2*rC2 - 16*alpha*re*rC -16*re2));  
    }


    public double getEpsilon() {return epsilon;}
 
    public final void setEpsilon(double eps) {
        epsilon = eps;
    }
    public Dimension getEpsilonDimension() {return Energy.DIMENSION;}
   

 
    public double getRe() {return re;}

    public final void setRe(double rEq) {
        re = rEq;
    }
    public Dimension getSigmaDimension() {return Length.DIMENSION;}
    
    
    
    public double getA() {return a;}

    public final void setA(double dummy) {
        a = dummy;
    }
    
    public double getZ1() {return z1;}

    public final void setZ1(double dummy) {
        z1 = dummy;
    }
    
    public double getZ2() {return z2;}

    public final void setZ2(double dummy) {
        z2 = dummy;
    }
    
    private static final long serialVersionUID = 1L;
    private double re;
    private double epsilon;
    private double a;
    private double z1;
    private double z2;
    private double coulombicCutoff;
}
