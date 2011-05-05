
package etomica.potential;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.BohrRadius;
import etomica.units.Hartree;

/**
 * Slim-lined version of P2HePCKLJS.java; main method compares values from each version.
 * 
 * Ab initio pair potential for helium from Przybytek et al. (2010) Phys. Rev. Lett. 104, 183003.  This is a true pair potential, rather than a pairwise-additive potential.
 * 
 * In this class, only the pair potential is valid, not the gradients, etc.  I am unlikely to ever include those...
 *
 * @author Kate Shaul
 */
public class P2HePCKLJSFancy extends Potential2SoftSpherical {
    
    public P2HePCKLJSFancy(ISpace space) {
        super(space);

        double W4 = 0.35322e-04/(alpha*alpha);
        double AS3= 0.577235e-06/(alpha*alpha*alpha);
        double polarizability = 1.38319217440;

        double K7 = 23.0/(4.0*Math.PI)*polarizability*polarizability/alpha;
        double ratio = alpha*K7/C6BO;

        B = new double[7];
        B[1] = 8.454943177941253;
        B[2] = 16.006586066260556;
        B[3] = 10.378373954734820;
        B[4] = 3.515803817223855;
        B[5] = 0.591502377533792;
        B[6] = 0.059455768329599;

        A = new double[7];
        A[1] = B[1];
        A[2] = B[2]-W4/C6BO;
        A[3] = B[3]-B[1]*W4/C6BO+AS3/C6BO;
        A[4] = ratio*B[5];
        A[5] = ratio*B[6];
        A[6] = 0;

        C = new double[17];
        C[3] = 0.000000577235; 
        C[4] = -0.000035322; 
        C[5] = 0.000001377841; 
        C[6] = 1.461830; 
        C[8] = 14.12350; 
        C[10] = 183.7497; 
        C[11] = -0.7674*100; 
        C[12] = 0.3372*1e4; 
        C[13] = -0.3806*1e4; 
        C[14] = 0.8534*1e5; 
        C[15] = -0.1707*1e6;  
        C[16] = 0.286*1e7; 
   
        P = new double[] { -25.4701669416621, 269.244425630616, -56.3879970402079};
        Q = new double[] {38.7957487310071, -2.76577136772754};
    }

    /**
     * The energy u.
     */
    public double u(double r2) {

    	double r = Math.sqrt(r2);
    	r = BohrRadius.UNIT.fromSim(r);

    	//Potential is speciously negative at separations less than 3 a0.
    	if (r < 0.3) {
    		return Double.POSITIVE_INFINITY;
    	}

    	double u1 = (P[0] + P[1]*r + P[2]*r*r)*Math.exp(-a*r);

    	double u2 = (Q[0] + Q[1]*r)*Math.exp(-b*r);

    	double br = eta*r;
        double m = Math.exp(-br);
        double term = m;
        sum[0] = term-1;
        for (int i=1; i<17; i++) {
            term *= br/i;
            sum[i] = sum[i-1] + term;
        }

    	double invr = 1.0/r;
    	double invr3 = invr*invr*invr;

    	double u3 = 0;
    	// jump in to the sum at i=3, invri = 1/r^3
    	double invri = invr3;
    	for (int i=3; i<17; i++) {
            u3 += sum[i]*C[i]*invri;
            invri *= invr;
    	}

    	/// damp_ret ////

        double sumA = 1.0;
        double sumB = 1.0;

        double x = alpha*r;
        double xn=1.0;
        for (int n=1;n<7;n++) {
            xn *= x;
            sumA += A[n]*xn;
            sumB += B[n]*xn;
        }
        double g = sumA/sumB;
    	double Vret = (C[3] + C[4]*invr + C6BO*(1.0-g)*invr3)*invr3;
        return Hartree.UNIT.toSim(u1 + u2 + u3 + Vret);
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
    	P2HePCKLJS p2Arrays = new P2HePCKLJS(space);
    	P2HePCKLJSFancy p2 = new P2HePCKLJSFancy(space);

    	//Test separations in fortran code:
    	double[] VTotRet = new double [] {286567.28801329771522433 ,    104344.43048737022036221 ,     36151.80990833066607593 ,     11964.29261818673148809 ,      3768.97799997742185951 ,      1111.08033532433000801 ,       292.72049836391101962 ,        58.44863402402420149 ,        -0.45682019466174945 ,        -9.15676266390127047 ,       -10.99557383520940412 ,        -9.67736339158196124 ,        -6.89071921436500023 ,        -4.62078614491897888 ,        -2.06561211026822855 ,        -0.98886027569716450 ,        -0.16554553396262922};
    	double[] rTest = new double [] {1.00 , 1.50 , 2.00 , 2.50 , 3.00 , 3.50 , 4.00 , 4.50 , 5.00 , 5.30 , 5.60 , 6.00 , 6.50 , 7.00 , 8.00 , 9.00 , 12.00};
    	double[] VTot = new double [] {286567.28368499479256570 ,    104344.42859618827060331 ,     36151.80886212403129321 ,     11964.29195944454113487 ,      3768.97754976155010809 ,      1111.08000968370220107 ,       292.72025283437358212 ,        58.44844291555266835 ,        -0.45697273095277996 ,        -9.15689724537757321 ,       -10.99569335167531925 ,        -9.67746633267005230 ,        -6.89080571621018390 ,        -4.62085971915516680 ,        -2.06566696229345581 ,        -0.98890251624988579 ,        -0.16556763818968243};
    	double[] VRet = new double [] {4.32830292764902966e-03 ,     1.89118193958213001e-03 ,     1.04620663571594000e-03 ,     6.58742191098119997e-04 ,     4.50215871792979995e-04 ,     3.25640627958039984e-04 ,     2.45529537389110018e-04 ,     1.91108471532230011e-04 ,     1.52536291030539990e-04 ,     1.34581476301949993e-04 ,     1.19516465915350001e-04 ,     1.02941088091070004e-04 ,     8.65018451842000018e-05 ,     7.35742361874499934e-05 ,     5.48520252273800021e-05 ,     4.22405527212700002e-05 ,     2.21042270531799994e-05};
    	double u; double uArrays;

    	System.out.println();
    	System.out.println("r(a0) \t V(Tot+Ret) no arrays \t V(Tot+Ret) arrays \t no arrays-arrays (K)   ");
     	for (int i=0;i<rTest.length;i++) {
    		double r = BohrRadius.UNIT.toSim(rTest[i]); //*AngstromPerBohrRadius; //Angstrom
    		u = p2.u(r*r); // Kelvin
    		uArrays = p2Arrays.u(r*r); // Kelvin

    		System.out.println(rTest[i]+"   \t"+u+"  \t" +uArrays+"  \t" +(u-uArrays) +" "+(u-uArrays)/(0.5*(u+uArrays)));
    	}

     	long t1 = System.currentTimeMillis();
        for (double r = 1; r<10; r+=0.0000001) {
            p2.u(r*r);
        }
        long t2 = System.currentTimeMillis();
        System.out.println((t2-t1)/1000.0);
    }

    private static final long serialVersionUID = 1L;
    protected final double[] sum = new double[17];
    protected final double[] C, A, B, P, Q;
    protected static final double C6BO = 1.460977837725; //pulled from potentials.f90...
    protected static final double alpha = 1.0/137.036; //fsalpha from potentials.f90...
    protected static final double a = 3.64890303652830;
    protected static final double b = 2.36824871743591;
    protected static final double eta = 4.09423805117871;
}
