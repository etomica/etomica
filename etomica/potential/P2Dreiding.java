package etomica.potential;

import etomica.space.Space;

/**
 * Dreiding potential.
 * 
 * U(r) = alpha * (R - Re)^2
 * 
 *  
 * where alpha is the potential's energy parameter, which equals 0.5*Ke (Ke = force constant)
 * 		 Re is the equilibrium bond radius
 * 
 * @author Tai Tan
 *
 */

public class P2Dreiding extends Potential2SoftSpherical {

	private double Re, alpha;
    
    public P2Dreiding(Space space, double Re, double alpha) {
        super(space);
        setAlpha(alpha);
        setBondLength(Re);
            }
 

    /**
     * harmonic oscillator
     */
    public double u(double r2) {
    	double r = Math.sqrt(r2);
    	double dr = r - Re;
        return alpha * dr * dr;
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
    	double r = Math.sqrt(r2);
    	//double diff = (r - Re)/Re*100;
    	//System.out.println(Re + " " + r + " " + diff);
        return 2*alpha*r*(r-Re);
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */

    public double d2u(double r2) {
        return 2*alpha*r2;
    }
            
    /**
     *  Integral used for corrections to potential truncation.
     */
    public double uInt(double rC) {
        return 0.0;
    }
    
    public double setAlpha(){return alpha;}
    public final void setAlpha(double a) {alpha = a;}
    
    
    public double setBondLength(){return Re;}
    public final void setBondLength (double b) {Re = b;}
 
	private static final long serialVersionUID = 1L; 
}
