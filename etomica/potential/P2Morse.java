package etomica.potential;
import etomica.EtomicaInfo;
import etomica.space.Space;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Length;

/**
 * Morse Potential interatomic potential
 * 
 *  The functional form of potential is :
 *  
 *  	u(r) = epsilon * {1- exp[-alpha * (r - re)]}^2			
 *   
 * where epsilon: describes the strength of the pair interaction
 *                (the well depth) 
 *         re   : is the equilibrium pair separation
 *        alpha : is a dimensionaless range parameter 
 *                (parameter controlling the width of the potential well)          
 *  
 *
 * @author Tai Boon Tan
 */
public final class P2Morse extends Potential2SoftSpherical {

    public P2Morse(Space space) {
        this(space, 1.0, 1.0, 1.0);
    }
    
    public P2Morse(Space space, double epsilon, double re, double alpha) {
        super(space);
        setEpsilon(epsilon);
        setRe(re);
        setAlpha(alpha);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Simple Lennard-Jones potential");
        return info;
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
    	double r = Math.sqrt(r2);
    	double expTerm1 = Math.exp(alpha*(1-(r/re)));
    	
    	return epsilon*(expTerm1-1)*(expTerm1-1)-epsilon;
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
    	double r = Math.sqrt(r2);
    	double expTerm = Math.exp(-alpha*(-1+(r/re)));
    	
    	return -2*(r/re)*epsilon*alpha*(expTerm-1)*expTerm;
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
    	double r = Math.sqrt(r2);
    	double expTerm = Math.exp(-alpha*(-1+(r/re)));
    	
        return 2*(r2/(re*re))*epsilon*alpha*alpha*expTerm*(2*expTerm-1);
    }
            
    /**
     *  Integral used for corrections to potential truncation.
     */
    public double uInt(double rC) {
        double A = space.sphereArea(1.0);  //multiplier for differential surface element
        double rC2 = rC*rC;
        double re2 = re*re;
        double alpha2 = alpha*alpha;
        double alpha3 = alpha*alpha*alpha;
        double expTerm = Math.exp(-alpha*(-1+(rC/re)));
        return (-A/(4*alpha3))
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
    
    
    
    public double getAlpha() {return alpha;}

    public final void setAlpha(double a) {
        alpha = a;
    }
    
    private static final long serialVersionUID = 1L;
    private double re;
    private double epsilon;
    private double alpha;
}
