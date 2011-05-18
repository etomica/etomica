
package etomica.potential;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.BohrRadius;
import etomica.units.Hartree;
import etomica.units.Kelvin;

/**
 * 
 * Ab initio pair potential for helium from Przybytek et al. (2010) Phys. Rev. Lett. 104, 183003.  This is a true pair potential, rather than a pairwise-additive potential.
 * 
 * In this class, only the pair potential is the second derivative is valid.  
 *
 * @author Kate Shaul
 */
public class P2HePCKLJS extends Potential2SoftSpherical {
    
    public P2HePCKLJS(ISpace space) {
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

    	//Potential is speciously negative at separations less than 0.3 a0.
    	if (r < 0.3) {
    		return Double.POSITIVE_INFINITY;
    	}

    	double u1 = (P[0] + P[1]*r + P[2]*r*r)*Math.exp(-a*r);

    	double u2 = (Q[0] + Q[1]*r)*Math.exp(-b*r);
    	
    	double invr = 1.0/r;
    	
    	double br = eta*r;
        double m = Math.exp(-br);
        double term = 1.0;
        double sum = term;
        double invri = invr;
        double u3 = 0;
        for (int i=1; i<17; i++) {
            term *= br/i;
            sum = sum + term;
            u3 += (-1.0+m*sum)*C[i]*invri;
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
        
        double invr3 = invr*invr*invr;
    	double Vret = (C[3] + C[4]*invr + C6BO*(1.0-g)*invr3)*invr3;
 
    	return  Hartree.UNIT.toSim(u1 + u2 + u3 + Vret);
    	
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        
    	double r = Math.sqrt(r2);
    	r = BohrRadius.UNIT.fromSim(r);

    	//Potential is speciously negative at separations less than 0.3 a0.
    	if (r < 0.3) {
    		return 0;
    	}

    	double u1 = (P[0] + P[1]*r + P[2]*r*r)*Math.exp(-a*r);
    	double du1dr = (P[1] + P[2]*2.0*r)*Math.exp(-a*r) - a*u1;
    	double u2 = (Q[0] + Q[1]*r)*Math.exp(-b*r);
    	double du2dr = Q[1]*Math.exp(-b*r) - b*u2;

    	
    	double invr = 1.0/r;

        double m = Math.exp(-eta*r);
        double dmdr = -eta*m;
        double term = 1.0;
        double sum = term;
        double dsumdr = 0;
        double u3 = 0;
    	double du3dr = 0;
        double invri = invr;
        for (int i=1; i<17; i++) {
            term *= eta*r/i;
            sum = sum + term;
            dsumdr = dsumdr + term*invr*i;
            u3 += (-1.0+m*sum)*C[i]*invri;
            du3dr += (dmdr*sum + m*dsumdr)*C[i]*invri + (-1.0+m*sum)*C[i]*invri*invr*(-i);
            invri *= invr;
        }


    	/// damp_ret ////

        double sumA = 1.0;
        double sumB = 1.0;
        double dsumAdr = 0;
        double dsumBdr = 0;
        double x = alpha*r;
        double xn=1.0;
        for (int n=1;n<7;n++) {
            xn *= x;
            sumA += A[n]*xn;
            sumB += B[n]*xn;
            dsumAdr += A[n]*xn*invr*n;
            dsumBdr += B[n]*xn*invr*n;
        }
        double g = sumA/sumB;
        double dgdr = dsumAdr/sumB - sumA/sumB/sumB*dsumBdr;
        double invr3 = invr*invr*invr;
    	//double Vret = (C[3] + C[4]*invr + C6BO*(1.0-g)*invr3)*invr3;
    	double dVretdr = (-3.0*C[3] -4.0*C[4]*invr -6.0*C6BO*(1.0-g)*invr3)*invr3*invr + C6BO*(-dgdr)*invr3*invr3;
        
    	double dudr = Hartree.UNIT.toSim(du1dr + du2dr + du3dr + dVretdr);
    	
        return r*dudr;
        
        
    
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
    	P2HePCKLJS p2 = new P2HePCKLJS(space);

    	//Test separations in fortran code, in Hartrees
    	double[] VTotRet = new double [] {9.07505678537836e-01,3.30439541259472e-01,1.14486105544985e-01,3.78887051832271e-02,1.19356572795740e-02,3.51858623016233e-03,9.26991759357220e-04,1.85096029792210e-04,-1.44666519197000e-06,-2.89977763063000e-05,-3.48209516983400e-05,-3.06464226675000e-05,-2.18216351894100e-05,-1.46331763646000e-05,-6.54141208063000e-06,-3.13153787265000e-06,-5.24252133480000e-07};
    	double[] rTest = new double [] {1.00 , 1.50 , 2.00 , 2.50 , 3.00 , 3.50 , 4.00 , 4.50 , 5.00 , 5.30 , 5.60 , 6.00 , 6.50 , 7.00 , 8.00 , 9.00 , 12.00}; 	
    	
    	double r;double u;
     	System.out.println();
    	System.out.println("r(a0) \t V(Tot+Ret) testEh.f90 \t V(Tot+Ret) here \t here-testEh.f90 (K)   ");
     	for (int i=0;i<rTest.length;i++) {
    		r = BohrRadius.UNIT.toSim(rTest[i]); //Angstrom
    		u = Hartree.UNIT.fromSim(p2.u(r*r)); // sim
    		double V = VTotRet[i];
    	 
    		System.out.println(rTest[i]+"   \t"+V+"  \t" +u+"  \t" +(u-V));
    	}

/*     	long t1 = System.currentTimeMillis();
        for (r = 1; r<10; r+=0.0000001) {
            p2.u(r*r);
        }
        long t2 = System.currentTimeMillis();
        System.out.println((t2-t1)/1000.0);*/
        
     	//Sanity check on rdu/dr
     	System.out.println();
     	r=5; 
     	double rmin = r;
     	double u1 = p2.u(r*r);
     	double delr = 0.000001;
     	System.out.println("r(A) \t\t rdudr \t\t\t rdeludelr \t\t\t\t u2");
     	while (r<(rmin+0.00001)) {

     		r = r + delr;
    		double u2 = p2.u(r*r); //u2 expects Angstroms - r is Angstroms
    		double rdudr = p2.du(r*r);
    		double rdeludelr = r*(u2-u1)/delr;
    		
    		//System.out.println(r+"   \t"+dudr+"  \t" +u2);
    		System.out.println(r+"   \t"+rdudr+"  \t" +rdeludelr+"  \t" +u2);
    		u1=u2;
    	}
    }
    
    private static final long serialVersionUID = 1L;
    protected final double[] C, A, B, P, Q;
    protected static final double C6BO = 1.460977837725; //pulled from potentials.f90...
    protected static final double alpha = 1.0/137.036; //fsalpha from potentials.f90...
    protected static final double a = 3.64890303652830;
    protected static final double b = 2.36824871743591;
    protected static final double eta = 4.09423805117871;

}
