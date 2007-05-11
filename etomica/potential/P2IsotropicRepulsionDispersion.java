package etomica.potential;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.simulation.Simulation;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.units.Dimension;
import etomica.units.Energy;

/**
 * Empirical Isotropic atom-atom repulsion-dispersion potential
 * Given formula:
 * 				
 * 				U(r) = A*exp(-r/B) - C /r^6	
 * 
 * 
 *	    A is in eV
 *      B is in Å
 *      C is in eV.Å^6
 *
 * @author Tai Tan
 */

public class P2IsotropicRepulsionDispersion extends Potential2SoftSpherical implements EtomicaElement {
	
	public P2IsotropicRepulsionDispersion(Simulation sim) {
        this(sim.getSpace(), sim.getDefaults().potentialWell, sim.getDefaults().atomSize, sim.getDefaults().potentialWell);
       
    }
	
    public P2IsotropicRepulsionDispersion(Space space, double AA, double BB, double CC) {
        super(space);
        dr01 = space.makeVector();
        setAA(AA);
        setBB(BB);
        setCC(CC);
        
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Dreiding Lennard-Jones potential");
        return info;
    }

    /**
     * The energy u
     */
    
    public double u(double r2) {
        double r = Math.sqrt(r2);
        
        if (r < 3*BB){
        	return Double.POSITIVE_INFINITY;
        }
        return AA*Math.exp(-r/BB) - CC/(r2*r2*r2);
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        double r = Math.sqrt(r2);
        return -AA/(BB)*r*Math.exp(-r/BB) + 6*CC/(r2*r2*r2);
        
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
        double r = Math.sqrt(r2);
        return AA/(BB*BB)*r2*Math.exp(-r/BB) - 42*CC/(r2*r2*r2);
    }
    
    /**
     *  Integral used for corrections to potential truncation.
     */
    public double uInt(double rC) { //need long range correction!!!!
        return 0.0;
    }

    public double getAA() {return AA;}
    public final void setAA(double a) {
        AA = a;
    }
    
    public double getBB() {return BB;}
    public final void setBB(double b) {
        BB = b;
        }
    
  
    public double getCC() {return CC;}
   public final void setCC(double c) {
        CC = c;
        }
   
    public Dimension getAADimension() {return Energy.DIMENSION;}
    public Dimension getBBDimension() {return Energy.DIMENSION;}
    public Dimension getCCDimension() {return Energy.DIMENSION;}
    
   
    private double AA, BB, CC;
    protected final IVector dr01;
	
	private static final long serialVersionUID = 1L;
}
