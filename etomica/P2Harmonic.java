package etomica;

/**
 * Harmonic Well interatomic potential.
 * Spherically symmetric potential of the form u(r) = 0.5*springConstant*(r-r0)^2 
 * where springConstant describes the strength of the pair interaction.
 *
 * @author Jhumpa Adhikari
 * @author David Kofke
 */

/* History 
 * 08/11/08 (DAK) added r0 parameter
 */
public class P2Harmonic extends Potential2SoftSpherical implements EtomicaElement {

    public String getVersion() {return "PotentialHarmonic:01.07.07/"+Potential2SoftSpherical.VERSION;}

    private double w = 100.0;// Spring constant gives a measure of the strength of harmonic interaction
	private final boolean r0Zero;
	private double r0;
    
    public P2Harmonic(double w) {
        this(Simulation.instance.hamiltonian.potential, w, 0.0);
    }
    public P2Harmonic(SimulationElement parent, double w) {
    	this(parent, w, 0.0);
    }
    /**
     * 
     * @param parent
     * @param w spring constant
     * @param r0  Separation at which potential is at its minimum.  Default is
     * zero.
     */
    public P2Harmonic(SimulationElement parent, double w, double r0) {
        super(parent);
        setSpringConstant(w);
        r0Zero = (r0 == 0.0);
        setR0(r0);
    }

    public double u(double r2) {
    	if(r0Zero) return 0.5*w*r2;
    	else {
    		double dr = Math.sqrt(r2) - r0;
    		return 0.5*w*dr*dr;
    	}
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
    	if(r0Zero) return w*r2;
    	else {
   			double r = Math.sqrt(r2);
   			return w*r*(r-r0);
    	}
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
    
	/**
	 * Separation at which potential is at its minimum.
	 * @return double
	 */
	public double getR0() {
		return r0;
	}

	/**
	 * Sets the the separation at which potential is at its minimum.
	 * @param r0 The r0 to set
	 */
	public void setR0(double r0) {
		this.r0 = r0;
	}

}//end of P2Harmonic
  