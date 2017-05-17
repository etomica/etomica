/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;
import etomica.space.Space;

/**
 * Pair potential for argon from Aziz (1993) JCP 99(6): 4518.  This is a true pair potential, rather than a pairwise-additive potential.
 * 
 * In this class, only the pair potential is valid, not the gradients, etc.  I am unlikely to ever include those...
 *
 * @author Kate Shaul
 */
public class P2ArgonAziz1993 extends Potential2SoftSpherical {
    
    public P2ArgonAziz1993(Space space) {
        super(space);
   
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
    	
    	double r = Math.sqrt(r2);
    	
    	double epsilonOverk = 143.235;  // Kelvin (convert to energy through Boltzmann's constant...)
    	double A = 87393.3927*epsilonOverk; // Kelvin
    	double alpha = 9.03228328/3.757; // inverse Angstroms
    	double beta = -0.168; // inverse Angstroms squared
        double uscf = A*Math.exp((-alpha*r)+(beta*r*r)); // Kelvin 
        
        r = r/(0.52917720859); // Bohr radius
        //r = r/(0.53);
        double rhor = 1.107*r;
        double f = 1.0 - Math.pow(rhor,1.68)*Math.exp(-0.78*rhor);
        double g6  =         63.5*Math.pow( (1.0-Math.exp(-2.1*rhor/6.0  - 0.109*rhor*rhor/Math.pow(6, 0.5) ))/r,  6); // Hartrees
        double g8  =       1510.0*Math.pow( (1.0-Math.exp(-2.1*rhor/8.0  - 0.109*rhor*rhor/Math.pow(8, 0.5) ))/r,  8); // Hartrees
        double g10 =      48000.0*Math.pow( (1.0-Math.exp(-2.1*rhor/10.0 - 0.109*rhor*rhor/Math.pow(10,0.5) ))/r, 10); // Hartrees
        double g12 =   2069581.26*Math.pow( (1.0-Math.exp(-2.1*rhor/12.0 - 0.109*rhor*rhor/Math.pow(12,0.5) ))/r, 12); // Hartrees
        double g14 =  116670633.0*Math.pow( (1.0-Math.exp(-2.1*rhor/14.0 - 0.109*rhor*rhor/Math.pow(14,0.5) ))/r, 14); // Hartrees
        double ucoor =-f*(g6 + g8 + g10 + g12 + g14)*(4.35974417E5/1.3806503); // Kelvin
        
        double u12 = uscf + ucoor;
        //System.out.println("");
        //System.out.println("reference u12/k (K) = " + u12);
        return uscf + ucoor; // Kelvin
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        
        return 0;
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
     
        return 0;
    }
            
    /**
     *  Integral used for corrections to potential truncation.
     */
    public double uInt(double rC) {
        
        return 0;  //complete LRC is obtained by multiplying by N1*N2/V
    }

   
   
    private static final long serialVersionUID = 1L;
    
}
