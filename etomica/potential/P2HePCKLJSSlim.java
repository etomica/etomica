
package etomica.potential;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Length;

/**
 * Slim-lined version of P2HePCKLJS.java; main method compares values from each version.
 * 
 * Ab initio pair potential for helium from Przybytek et al. (2010) Phys. Rev. Lett. 104, 183003.  This is a true pair potential, rather than a pairwise-additive potential.
 * 
 * In this class, only the pair potential is valid, not the gradients, etc.  I am unlikely to ever include those...
 *
 * @author Kate Shaul
 */
public class P2HePCKLJSSlim extends Potential2SoftSpherical {
    
    public P2HePCKLJSSlim(ISpace space) {
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

    	double u1 = (-25.4701669416621+ 269.244425630616*r+ -56.3879970402079*r*r)*Math.exp(-3.64890303652830*r);
    	
    	double u2 = (38.7957487310071 + -2.76577136772754*r)*Math.exp(-2.36824871743591*r);

    	double br = 4.09423805117871*r;	
    	double term = br*br/2.0;
    	double sum2 = 1.0 + br + term;
    	term = term*br/3.0;
    	double sum3 = sum2 + term;
    	term = term*br/4.0;
    	double sum4 = sum3 + term;
    	term = term*br/5.0;
    	double sum5 = sum4 + term;
    	term = term*br/6.0;
    	double sum6 = sum5 + term;
    	term = term*br/7.0;
    	double sum7 = sum6 + term;
    	term = term*br/8.0;
    	double sum8 = sum7 + term;
    	term = term*br/9.0;
    	double sum9 = sum8 + term;
    	term = term*br/10.0;
    	double sum10 = sum9 + term;
    	term = term*br/11.0;
    	double sum11 = sum10 + term;
    	term = term*br/12.0;
    	double sum12 = sum11 + term;
    	term = term*br/13.0;
    	double sum13 = sum12 + term;
    	term = term*br/14.0;
    	double sum14 = sum13 + term;
    	term = term*br/15.0;
    	double sum15 = sum14 + term;
    	term = term*br/16.0;
    	double sum16 = sum15 + term;
    	
    	double invr2 = 1.0/(r*r);
    	double invr3 = invr2/r;
    	double invr4 = invr2*invr2;
    	double invr5 = invr3*invr2;
    	double invr6 = invr3*invr3;
    	double invr7 = invr3*invr4;
    	double invr8 = invr4*invr4;
    	double m = Math.exp(-br);

    	double u3 = -(1.0-m*sum3 )*0.000000577235*invr3   -(1.0-m*sum4)*-0.000035322*invr4    -(1.0-m*sum5)*0.000001377841*invr5;
    	u3 = u3     -(1.0-m*sum6 )*1.461830*invr6   -(1.0-m*sum8)*14.12350*invr8   -(1.0-m*sum10)*183.7497*invr5*invr5;
    	u3 = u3     -(1.0-m*sum11)*-0.7674*100*invr5*invr6 -(1.0-m*sum12)*0.3372*1e4*invr6*invr6 -(1.0-m*sum13)*-0.3806*1e4*invr5*invr8 ;
    	u3 = u3     -(1.0-m*sum14)*0.8534*1e5*invr7*invr7 -(1.0-m*sum15)*-0.1707*1e6*invr7*invr8 -(1.0-m*sum16)*0.286*1e7*invr8*invr8;

    	
    	/// damp_ret ////

        double x = r/137.036;        

        double g = (1.0+x*(8.454943177941253+x*(15.552570226418524+x*(7.5564431748600205+x*(1.4177376898763494+x*0.1425060774782998)))))/(1.0+x*(8.454943177941253+x*(16.006586066260556+x*(10.378373954734820+x*(3.515803817223855+x*(0.591502377533792+x*0.059455768329599))))));
        //double g = (1.0+x*(A1+x*(A2+x*(A3+x*(A4+x*A5)))))/(1.0+x*(B1+x*(B2+x*(B3+x*(B4+x*(B5+x*B6))))));
    	double Vret = 0.000000577235*invr3 + -0.000035322*invr4 + 1.460977837725*(1.0-g)*invr6;
       
        return  (u1 + u2 + u3 + Vret)*KPerHartree; // Kelvin

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
    	P2HePCKLJSSlim p2 = new P2HePCKLJSSlim(space);
    	
    	//Test separations in fortran code:
    	double[] VTotRet = new double [] {286567.28801329771522433 ,    104344.43048737022036221 ,     36151.80990833066607593 ,     11964.29261818673148809 ,      3768.97799997742185951 ,      1111.08033532433000801 ,       292.72049836391101962 ,        58.44863402402420149 ,        -0.45682019466174945 ,        -9.15676266390127047 ,       -10.99557383520940412 ,        -9.67736339158196124 ,        -6.89071921436500023 ,        -4.62078614491897888 ,        -2.06561211026822855 ,        -0.98886027569716450 ,        -0.16554553396262922};
    	double[] rTest = new double [] {1.00 , 1.50 , 2.00 , 2.50 , 3.00 , 3.50 , 4.00 , 4.50 , 5.00 , 5.30 , 5.60 , 6.00 , 6.50 , 7.00 , 8.00 , 9.00 , 12.00};
    	double[] VTot = new double [] {286567.28368499479256570 ,    104344.42859618827060331 ,     36151.80886212403129321 ,     11964.29195944454113487 ,      3768.97754976155010809 ,      1111.08000968370220107 ,       292.72025283437358212 ,        58.44844291555266835 ,        -0.45697273095277996 ,        -9.15689724537757321 ,       -10.99569335167531925 ,        -9.67746633267005230 ,        -6.89080571621018390 ,        -4.62085971915516680 ,        -2.06566696229345581 ,        -0.98890251624988579 ,        -0.16556763818968243};
    	double[] VRet = new double [] {4.32830292764902966e-03 ,     1.89118193958213001e-03 ,     1.04620663571594000e-03 ,     6.58742191098119997e-04 ,     4.50215871792979995e-04 ,     3.25640627958039984e-04 ,     2.45529537389110018e-04 ,     1.91108471532230011e-04 ,     1.52536291030539990e-04 ,     1.34581476301949993e-04 ,     1.19516465915350001e-04 ,     1.02941088091070004e-04 ,     8.65018451842000018e-05 ,     7.35742361874499934e-05 ,     5.48520252273800021e-05 ,     4.22405527212700002e-05 ,     2.21042270531799994e-05};
    	double u; double uArrays;
    	

    	System.out.println();
    	System.out.println("r(a0) \t V(Tot+Ret) no arrays \t V(Tot+Ret) arrays \t no arrays-arrays (K)   ");
     	for (int i=0;i<rTest.length;i++) {
    		double r = rTest[i]*AngstromPerBohrRadius; //Angstrom
    		u = p2.u(r*r); // Kelvin
    		uArrays = p2Arrays.u(r*r); // Kelvin
    		
    	 
    		System.out.println(rTest[i]+"   \t"+u+"  \t" +uArrays+"  \t" +(u-uArrays));
    	}
     	

    }


    private static final double AngstromPerBohrRadius = 0.529177; // Rounding provided by Pryzbytek et al. 2010
    private static final double KPerHartree = 315774.65; // Rounding provided by Pryzbytek et al. 2010
    private static final long serialVersionUID = 1L;
    
}
