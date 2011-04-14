
package etomica.potential;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Length;

/**
 * Ab initio pair potential for helium from Przybytek et al. (2010) Phys. Rev. Lett. 104, 183003.  This is a true pair potential, rather than a pairwise-additive potential.
 * 
 * In this class, only the pair potential is valid, not the gradients, etc.  I am unlikely to ever include those...
 *
 * @author Kate Shaul
 */
public class P2HePCKLJ extends Potential2SoftSpherical {
    
    public P2HePCKLJ(ISpace space) {
        super(space);
   
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
    	
    	double r = Math.sqrt(r2);

    	r = r/AngstromPerBohrRadius; // Bohr radius
    	
    	//Potential is speciously negative at separations less than 3 a0.
    	if (r < 0.3) {
    		return Double.POSITIVE_INFINITY;
    	}
    	
    	//Parameters
    	
    	double[] C = new double[17];
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
    	
    	double C6BO = 1.460977837725; //pulled from potentials.f90...
    	
    	double[] APaper = new double[6];
    	APaper[1] = 8.454943177941253;
    	APaper[2] = 15.552891567112597; 
    	APaper[3] = 7.559160092169168; 
    	APaper[4] = 1.417737689876350; 
    	APaper[5] = 0.142506077478301; 
    	
    	

    	
    	double[] B = new double[7];
    	B[1] = 8.454943177941253;
    	B[2] = 16.006586066260556;
    	B[3] = 10.378373954734820;
    	B[4] = 3.515803817223855;
    	B[5] = 0.591502377533792;
    	B[6] = 0.059455768329599;

    	double alpha = 1.0/137.036; //fsalpha from potentials.f90...
    	double W4 = 0.35322e-04/(alpha*alpha);
    	double AS3= 0.577235e-06/(alpha*alpha*alpha);
    	double polarizability = 1.38319217440;

    	double K7 = 23.0/(4.0*Math.PI)*polarizability*polarizability/alpha;
    	double ratio = alpha*K7/C6BO;

    	double[] A = new double[6];
    	A[1] = B[1];
        A[2] = B[2]-W4/C6BO;
        A[3] = B[3]-B[1]*W4/C6BO+AS3/C6BO;
        A[4] = ratio*B[5];
        A[5] = ratio*B[6];

    	double a = 3.64890303652830;
    	double b = 2.36824871743591;
    	double eta = 4.09423805117871;
    	
    	double[] P = new double[] { -25.4701669416621, 269.244425630616, -56.3879970402079};
    	double[] Q = new double[] {38.7957487310071, -2.76577136772754};

    	

    	double u1 = 0;
    	for (int i=0; i<=2; i++) {
    		
    		u1 = u1 + P[i]*Math.pow(r,i);
    		
    	}
    	u1 = u1*Math.exp(-a*r);
    	
    	double u2 = 0;
    	for (int i=0; i<=1; i++) {
    		
    		u2 = u2 + Q[i]*Math.pow(r,i);
    		
    	}
    	u2 = u2*Math.exp(-b*r);
    	
    	/// damp_ret ////
        double sumA = 0;
        double sumB = 0;
        
        double x = alpha*r;
        double xn=1.0;
        for (int n=1;n<=5;n++) {
        	xn=xn*x;
        	sumA = sumA + A[n]*xn;
        	sumB = sumB + B[n]*xn;
        } 
        sumB = sumB + B[6]*xn*x;
        double g = (1.0+sumA)/(1.0+sumB);
        
    	double u3 = 0;
    	double Vret = 0;
    	double br = eta*r;
    	double sum = 1.0 + br + br*br/2.0;
    	double term = br*br/2.0;
    	double f = 1.0;
    	double asy = 0;
    	double asy2 = 0;
    	for (int n=3; n<=16; n++) {
    		
    		term = term*br/n;
    		sum = sum + term;
    		
    		//issue of br>1.0 or <1.0 not in paper.
    		if (br>1.0) {
    			f = 1.0-(Math.exp(-br)*sum); 
    		} else {
    			
    			double suma=0;
    			double terma = term;
    			for (int i=n+1; i<=n+40; i++) {
    				terma = terma*br/n;
    				suma=suma+term;
    			}
    			f = suma*Math.exp(-br); 
    		}
    		
 		
    		u3 = u3 - f*C[n]/Math.pow(r,n);
    	    
    		double fmod = -Math.exp(-br)*sum;
    		
    	    if (n<5) {  	    	
    	    	Vret = Vret + f*C[n]/Math.pow(r,n);
    	    	Vret = Vret - fmod*C[n]/Math.pow(r,n);
    	    	
    	    	
    	    	
    	    } 
    	    
    	    if (n==6) {
    	    	Vret = Vret + f*C[n]/Math.pow(r,6);

    	    	Vret = Vret - (g+fmod)*1.460977837725/Math.pow(r,6);
    	    	Vret = Vret - f*(C[n]-1.460977837725)/Math.pow(r,6);

    	    	
    	    } 
    		
    	}

    	
 
        double u = u1 + u2 + u3; // Hartree
       
        //System.out.println(asy +" "+ asy2+" "+ (asy+asy2)*KPerHartree);
        //System.out.println(g);
        if (part == 1) {
        	return u*KPerHartree; // Kelvin
        } else if (part == 2) {
        	return Vret*KPerHartree; // Kelvin
        } else if (part == 3) {
        	return  (u+Vret)*KPerHartree; // Kelvin
        } else {
        	throw new RuntimeException("Part can be only 1,2,or3");
        }

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
    	P2HePCKLJ p2 = new P2HePCKLJ(space);
    	
    	//Test separations in fortran code:
    	double[] VTotRet = new double [] {286567.28801329771522433 ,    104344.43048737022036221 ,     36151.80990833066607593 ,     11964.29261818673148809 ,      3768.97799997742185951 ,      1111.08033532433000801 ,       292.72049836391101962 ,        58.44863402402420149 ,        -0.45682019466174945 ,        -9.15676266390127047 ,       -10.99557383520940412 ,        -9.67736339158196124 ,        -6.89071921436500023 ,        -4.62078614491897888 ,        -2.06561211026822855 ,        -0.98886027569716450 ,        -0.16554553396262922};
    	double[] rTest = new double [] {1.00 , 1.50 , 2.00 , 2.50 , 3.00 , 3.50 , 4.00 , 4.50 , 5.00 , 5.30 , 5.60 , 6.00 , 6.50 , 7.00 , 8.00 , 9.00 , 12.00};
    	double[] VTot = new double [] {286567.28368499479256570 ,    104344.42859618827060331 ,     36151.80886212403129321 ,     11964.29195944454113487 ,      3768.97754976155010809 ,      1111.08000968370220107 ,       292.72025283437358212 ,        58.44844291555266835 ,        -0.45697273095277996 ,        -9.15689724537757321 ,       -10.99569335167531925 ,        -9.67746633267005230 ,        -6.89080571621018390 ,        -4.62085971915516680 ,        -2.06566696229345581 ,        -0.98890251624988579 ,        -0.16556763818968243};
    	double[] VRet = new double [] {4.32830292764902966e-03 ,     1.89118193958213001e-03 ,     1.04620663571594000e-03 ,     6.58742191098119997e-04 ,     4.50215871792979995e-04 ,     3.25640627958039984e-04 ,     2.45529537389110018e-04 ,     1.91108471532230011e-04 ,     1.52536291030539990e-04 ,     1.34581476301949993e-04 ,     1.19516465915350001e-04 ,     1.02941088091070004e-04 ,     8.65018451842000018e-05 ,     7.35742361874499934e-05 ,     5.48520252273800021e-05 ,     4.22405527212700002e-05 ,     2.21042270531799994e-05};
    	double u;
    	
    	p2.setPart(3);
    	System.out.println();
    	System.out.println("r(a0) \t V(Tot+Ret) test.f90 \t V(Tot+Ret) here \t here-test.f90 (K)   ");
     	for (int i=0;i<rTest.length;i++) {
    		double r = rTest[i]*AngstromPerBohrRadius; //Angstrom
    		u = p2.u(r*r); // Kelvin
    		
    	 
    		System.out.println(rTest[i]+"   \t"+VTotRet[i]+"  \t" +u+"  \t" +(u-VTotRet[i]));
    	}
     	
     	p2.setPart(1);
     	System.out.println();
    	System.out.println();
    	System.out.println("r(a0) \t VTot (K) test.f90 \t VTot (K) here  \t here-test.f90 (K)  ");
     	for (int i=0;i<rTest.length;i++) {
    		double r = rTest[i]*AngstromPerBohrRadius; //Angstrom
    		
    		u = p2.u(r*r); // Kelvin
    		
    	 
    		System.out.println(rTest[i]+"   \t"+VTot[i]+"  \t" +u+"  \t" +(u-VTot[i]));
    	}
     	
     	p2.setPart(2);
     	System.out.println();
    	System.out.println();
    	System.out.println("r(a0) \t VRet (K) (test.f90) \t VRet (K) (here) \t here-test.f90 (K)   ");
     	for (int i=0;i<rTest.length;i++) {
    		double r = rTest[i]*AngstromPerBohrRadius; //Angstrom
    		u = p2.u(r*r); // Kelvin
    		
    	 
    		System.out.println(rTest[i]+"   \t"+VRet[i]+"  \t" +u + "  \t" + (u-VRet[i]));
    	}

    	
    }
   
    public void setPart(int i) {
    	part = i;
    }

    private static final double AngstromPerBohrRadius = 0.529177; // Rounding provided by Pryzbytek et al. 2010
    private static final double KPerHartree = 315774.65; // Rounding provided by Pryzbytek et al. 2010
    private static final long serialVersionUID = 1L;
    private int part = 3; // default is 3: both Vtot(1) and Vret(2)
    
}
