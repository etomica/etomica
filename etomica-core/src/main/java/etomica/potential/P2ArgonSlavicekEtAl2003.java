/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;
import etomica.space.Space;
import etomica.space3d.Space3D;

/**
 * Ab initio pair potential for argon from Slavicek et al. (1993) JCP 119(4): 2102.  This is a true pair potential, rather than a pairwise-additive potential.
 * <p>
 * In this class, only the pair potential is valid, not the gradients, etc.  I am unlikely to ever include those...
 *
 * @author Kate Shaul
 */
public class P2ArgonSlavicekEtAl2003 extends Potential2SoftSpherical {

	public static Potential2Soft makeTruncated(Space space, TruncationFactory tf) {
		return tf.make(new P2ArgonSlavicekEtAl2003(space));
	}

	public P2ArgonSlavicekEtAl2003(Space space) {
		super(space);

	}

	/**
	 * The energy u.
	 */
	public double u(double r2) {
    	
    	double r = Math.sqrt(r2);
    
    	double r6 = r2*r2*r2;
    	double r8 = r6*r2;
    	double r10 = r8*r2;
    	
    	
    	double rmin = 3.771; //Angstroms
    	double rmin2 = rmin*rmin;
    	double rmin6 = rmin2*rmin2*rmin2;
    	double rmin8 = rmin6*rmin2;
    	double rmin10 = rmin8*rmin2;
    	
    	double eOverkB = 142.331; // Kelvin
    	
    	//Repulsive (SCF) component
    	
    	double A = 333634.40*eOverkB; 
    	double a = 11.493139/rmin;
    	double b = -1.082649/rmin2;
    	
    	double uSCF = A*Math.exp(-a*r + b*r*r);
    		
    	
    	// Dispersion (correlation) component
    	double C6 = 1.1261753*eOverkB*rmin6;
    	double C8 = 0.34453349*eOverkB*rmin8;
    	double C10 = 0.72590255*eOverkB*rmin10;
    	double D = 1.1422763*rmin;
    	
    	double F = 1.0;
    	
    	if (r < D) {
    		F = Math.exp(-((D/r) -1)*((D/r) -1) );
    	}
    	
    	double uCOR = -F*( (C6/r6)+(C8/r8)+(C10/r10) );
    	
        double u = uSCF+uCOR;// Kelvin
       
        return u; // Kelvin
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
    
    public static void main(String[] args) {
    	
    	Space space = Space3D.getInstance();
    	P2ArgonSlavicekEtAl2003 p2 = new P2ArgonSlavicekEtAl2003(space);
    	
    	//Minimum of the potential published by Slavicek et al. (2003):
   
    	double r = 3.6;
    	double u;
    	
    	
    	double umin = 0;
    	double rmin = 0;
    	while (r<4) {
    		r = r + 0.001;
    		u = p2.u(r*r); // Kelvin
    		u = u/KPerHartree*1e6;// microHartree
    		if (u < umin) {
    			
    			umin = u;
    			rmin =r;
    		}
    		//System.out.println(r+"  "+u);
    	}
    	
    	//Minimum
    	System.out.println(rmin+"  "+umin);
    	
    	
    	r = 3.757178;
		u = p2.u(r*r); // Kelvin
		u = u/KPerHartree*1e6;// microHartree
		System.out.println(r+"  "+u);
		
    	
    }
   

    private static final double rBohr = 0.5292; // Rounding provided by Tang and Toennies
    private static final double KPerHartree = 3.158e5; // Rounding provided by Tang and Toennies
    //private static final double kB =  1.3806503e-23;// 1.38e-23; // 1.3806503e-23 J/Kelvin
    //private static final double JPerHartree = 4.359743e-18; // 4.359743e-18 J;
    //private static final double KPerHartree = 4.359743e-18/1.3806503e-23; // Rounding provided by Tang and Toennies 
    private static final long serialVersionUID = 1L;
    
}
