package etomica.potential;

import etomica.EtomicaElement;
import etomica.Space;


/**
 * Harmonic Well interatomic potential.
 * Spherically symmetric potential of the form u(r) = 0.5*springConstant*(r)^2
 * where springConstant describes the strength of the pair interaction.
 *
 * @author Jhumpa Adhikari
 * @author David Kofke
 */
public class P2Anharmonic extends Potential2SoftSpherical implements EtomicaElement {

    private double w = 100.0;// Spring constant gives a measure of the strength of harmonic interaction
    private double a = 0;	// Anharmonic spring constant

    public P2Anharmonic(Space space, double w) {
        super(space);
        setSpringConstant(w);
    }

    public double u(double r2) {
        return r2*(0.5*w+a*Math.sqrt(r2));
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        return r2*(w+3*a*Math.sqrt(r2));
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
        return r2*(w+6*a*Math.sqrt(r2));
    }
            
    /**
     *  Integral used for corrections to potential truncation.
     */
    public double uInt(double rC) {
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
    public double getAnharmonicConstant() {return a;}
    public void setAnharmonicConstant(double factor) {a=factor;}
}//end of P2Harmonic
  