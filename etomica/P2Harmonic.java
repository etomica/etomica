package etomica;

/**
 * Harmonic Well interatomic potential.
 * Spherically symmetric potential of the form u(r) = 0.5*springConstant*(r)^2
 * where springConstant describes the strength of the pair interaction.
 *
 * @author Jhumpa Adhikari
 * @author David Kofke
 */
public class P2Harmonic extends Potential2SoftSpherical implements EtomicaElement {

    public String getVersion() {return "PotentialHarmonic:01.07.07/"+Potential2SoftSpherical.VERSION;}

    private double w = 100.0;// Spring constant gives a measure of the strength of harmonic interaction
    
    public P2Harmonic(double w) {
        this(Simulation.instance, w);
    }
    public P2Harmonic(Simulation sim, double w) {
        super(sim);
        setSpringConstant(w);
    }

    public double u(double r2) {
        return 0.5*w*r2;
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        return w*r2;
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
        return w*r2;
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
    
}//end of P2Harmonic
  