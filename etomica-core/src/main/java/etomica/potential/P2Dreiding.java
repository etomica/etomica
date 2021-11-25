/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

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

public class P2Dreiding implements Potential2Soft {

    public static Potential2Soft makeTruncated(double Re, double alpha, TruncationFactory tf) {
        return tf.make(new P2Dreiding(Re, alpha));
    }

    private double Re, alpha;

    public P2Dreiding(double Re, double alpha) {
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
        return 2*alpha*r*(r-Re);
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */

    public double d2u(double r2) {
        return 2*alpha*r2;
    }

    public double getAlpha(){return alpha;}
    public final void setAlpha(double a) {alpha = a;}
    
    
    public double getBondLength(){return Re;}
    public final void setBondLength (double b) {Re = b;}
 
	private static final long serialVersionUID = 1L; 
}
