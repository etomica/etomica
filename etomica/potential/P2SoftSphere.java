package etomica.potential;
import etomica.EtomicaInfo;
import etomica.space.ISpace;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Length;

/**
 * Soft-sphere interatomic potential.
 * Characterized by the pairwise-additive, spherically symmetric,
 *  inverse-power potential of the form:
 *  
 * 			     u(r) = epsilon*[sigma/ r]^n 
 * 
 * where epsilon:  describes the strength of the pair interaction, 
 *       sigma  :  is the atom size parameter
 *         n    :  is the degree of softness, s=1/n, 
 *               eg: hard sphere when s=0 or n->infinity
 *               
 *     sigma/r  : sig_r
 * 
 *
 * @author Tai Boon Tan
 */
public final class P2SoftSphere extends Potential2SoftSpherical {

    public P2SoftSphere(ISpace space) {
        this(space, 1.0, 1.0, 12);

    }
    
    public P2SoftSphere(ISpace space, double sigma, double epsilon, int n) {

        super(space);
        setSigma(sigma);
        setEpsilon(epsilon);
        setExponent(n);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Soft-Sphere potential");
        return info;
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
    	double sig_r = sigma/Math.sqrt(r2);
    	int expN = getExponent();
    	double sig_rn = 1;
    	
    	for (int i=0; i<expN; i++){
    			sig_rn *= sig_r;
    	}
    	
    	return epsilon*sig_rn;
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
    	double sig_r = sigma/Math.sqrt(r2);
    	int expN = getExponent();
    	double sig_rn = 1;
    
    	for (int i=0; i<expN; i++){
    				sig_rn *= sig_r;
    	}
    	
        return -n*epsilon*sig_rn;
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
    	double sig_r = sigma/Math.sqrt(r2);
    	int expN = getExponent();
    	double sig_rn = 1;
    
    	for (int i=0; i<expN; i++){
    				sig_rn *= sig_r;
    	}
        return n*(n+1)*epsilon*sig_rn;
    }
            
    /**
     *  Integral used for corrections to potential truncation.
     */
    public double uInt(double rC) {
        double A = space.sphereArea(1.0);  //multiplier for differential surface element
        int D = space.D();                 //spatial dimension
        double rC3 = rC*rC*rC;
         
    	int expN = getExponent();
    	double sig_rCn = sigma/rC;
    	
    	if (expN ==0){
    		sig_rCn = 1;
    	} else
    		{
    			for (int i=0; i<expN; i++){
    				sig_rCn *= sigma/rC;
    			}
    		}
       return epsilon*A*rC3*sig_rCn/(n-D);  //complete LRC is obtained by multiplying by N1*N2/V
    }

    /**
     * Accessor method for soft-sphere size parameter.
     */
    public double getSigma() {return sigma;}
    /**
     * Mutator method for soft sphere size parameter.
     * Does not adjust potential cutoff for change in sigma.
     */
    public final void setSigma(double sig) {
        sigma = sig;
    }
    public Dimension getSigmaDimension() {return Length.DIMENSION;}
    
    /**
     * Accessor method for soft-sphere energy parameter
     */
    public double getEpsilon() {return epsilon;}
    /**
     * Mutator method for soft-sphere energy parameter
     */
    public final void setEpsilon(double eps) {
        epsilon = eps;
    }
    public Dimension getEpsilonDimension() {return Energy.DIMENSION;}
   
    /**
     * Accessor method for soft-sphere softness parameter
     */
    public int getExponent() {return n;}
    /**
     * Mutator method for soft-sphere softness parameter
     */
    public final void setExponent(int en) {
        n = en;
    }
    
    private static final long serialVersionUID = 1L;
    private double sigma;
    private double epsilon;
    private int n;
}
