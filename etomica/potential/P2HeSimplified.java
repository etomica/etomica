
package etomica.potential;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.space3d.Space3D;

/**
 * Simplified pair potential for Helium.  The potential has the form of an
 * exponential-6 model with an additional r^-10 term.
 *
 * @author Andrew Schultz
 */
public class P2HeSimplified extends Potential2SoftSpherical {
    
    public P2HeSimplified(ISpace space) {
        super(space);
    }

    /**
     * The energy u.
     */
    public double u(double r2) {

    	if (r2 < sigmaHC2) {
    		return Double.POSITIVE_INFINITY;
    	}

    	double r = Math.sqrt(r2);
        double r4 = r2*r2;
        double r6 = r4*r2;
        double r10 = r6*r4;
        return A0*Math.exp(-A1*r)-A2/r6-A3/r10;
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
        return 0;
    }
    
    
    
    public static void main(String[] args) {
        // this main method simply performs a benchmark for this class
        // against P2HePCKLJS (the full Helium potential).
    	Space space = Space3D.getInstance();
    	P2HeSimplified p2Simple = new P2HeSimplified(space);
    	P2HePCKLJS p2 = new P2HePCKLJS(space);

    	int n = 10000000;
    	long t1=System.currentTimeMillis();
    	double totU = 0;
        double rr2 = 9;
    	for (int irr=0; irr<n; irr++) {
    	    p2.u(rr2);
    	}
        for (int irr=0; irr<n; irr++) {
            p2.u(rr2);
        }
        System.out.println(totU);
        long t2 = System.currentTimeMillis();

        long t1s=System.currentTimeMillis();
        for (int irr=0; irr<n; irr++) {
            p2Simple.u(rr2);
        }
        for (int irr=0; irr<n; irr++) {
            p2Simple.u(rr2);
        }
        long t2s = System.currentTimeMillis();
        System.out.println("PCKLJS "+(t2-t1)/1000.0+"     simple "+(t2s-t1s)/1000.0+"     "+totU);
    }
    
    private static final long serialVersionUID = 1L;
    protected final double A0 = 9.05469e+06;
    protected final double A1 = 4.60617;
    protected final double A2 = 9021.61;
    protected final double A3 = 339994;
    protected final double sigmaHC2 = 1.7*1.7;
}
