package etomica.math;

/**
 * Static-method library of various functions
 */
 
public final class SpecialFunctions {
    
    private SpecialFunctions() {}
    
    
    /**
     * Complementary error function, computed using the approximant 7.1.26 of Abramowitz & Stegun.
     * Defined for x >= 0
     */
    public static double erfc(double x) {
        double t = 1.0/(1.0 + 0.3275911*x);
        return (t * (
                  0.254829592 + t * (
                     -0.284496736 + t * (
                       1.421413741 +  t * (
                         - 1.453152027 + 1.061405429*t)))))*Math.exp(-x*x);
    }
    
    /**
     * The factorial function, n!
     * 
     * @throws IllegalArgumentException if n < 0
     */
    public static int factorial(int n){
        if(n < 0){
            throw new IllegalArgumentException("Argument less than zero: "+n);
        }
        return (n <= 1) ? 1 :(n*factorial(n-1));
    }
    
    public static double lnFactorial(int n) {
        if(n < 0) {
            throw new IllegalArgumentException("Argument less than zero: "+n);
        }
        return (n <= 1) ? 0 :(Math.log(n) + lnFactorial(n-1));
    }
        
    //non-recursive version
//	public static int factorial(int n) {
//		if(n < 0) throw new IllegalArgumentException("Illegal to pass negative value to factorial");
//		int factorial = 1;
//		for (int i=2; i<=n; i++) {
//		   factorial *= i;
//		}
//		return factorial;	
//	}

    /**
     * The sign function, returning -1 if x < 0, zero if x == 0, and +1 if x > 0.
     */
    public static double sgn(double x) {
        return (x < 0.0) ? -1.0 : ((x > 0.0) ? +1.0 : 0.0);
    }


    /**
     * Returns the ln(gamma), the natural logarithm of the gamma function.
     * This method is not tested.
     */

    public static double lnGamma(double x) {
    	double tmp = x+5.5;
    	double y = x;
    	double stp = 2.5066282746310005;
    	tmp = (x+0.5)*Math.log(tmp) - tmp;
    	double ser = 1.000000000190015;
    	for(int i=0; i<5; i++) {
    		y += 1.0;
    		ser += lnGammaCoeff[i]/y;
    	}
    	return tmp + Math.log(stp*ser/x);
    }
    private static final double[] lnGammaCoeff = new double[] {
    		+76.18009172947146,
			-86.50532032941677,
			+24.01409824083091,
			-1.231739572450155,
			+0.1208650973866179e-02,
			-0.5395239384953e-05};
    
    /**
     * Normalized incomplete gamma function, equal to
     * Integrate[Exp[-t] t^(a-1),{t,x,Infinity}]/Gamma[a]
     */
    //method is not tested
    public static double gammaQ(double a, double x) {
    	if(x < 0.0 || a < 0.0) throw new IllegalArgumentException();
    	if(a == 0.0) return 1.0;
    	if(x < a+1.0) {
    		return 1.0 - gser(a, x);
    	} else {
    		return gcf(a, x);
    	}
    }
    
    private static double gser(double a, double x) {
    	int nmax = 500;
    	double epsilon = 3.0e-12;
    	double ap = a;
    	double sum = 1.0/a;
    	double del = sum;
    	for(int n=1; n<=nmax; n++) {
    		ap++;
    		del *= x/ap;
    		sum += del;
    		if(Math.abs(del) < Math.abs(sum)*epsilon) break;
    	}
    	return sum * Math.exp(-x + a*Math.log(x) - lnGamma(a));
    }
    
    private static double gcf(double a, double x) {
    	int imax = 500;
    	double epsilon = 3.0e-12;
    	double b = x + 1.0 - a;
    	double c = 1.0/Double.MIN_VALUE;
    	double d = 1.0/b;
		double h = d;
    	for(int i=0; i<=imax; i++) {
    		double an = -i*(i-a);
    		b += 2.0;
    		d = an*d + b;
    		if(Math.abs(d) < Double.MIN_VALUE) d = Double.MIN_VALUE;
    		c = b + an/c;
    		if(Math.abs(c) < Double.MIN_VALUE) c = Double.MIN_VALUE;
    		d = 1.0/d;
    		double del = d*c;
    		h *= del;
    		if(Math.abs(del-1.0) < epsilon) break;
    	}
    	return Math.exp(-x + a*Math.log(x) - lnGamma(a))*h;
    }
    
    public static void main(String[] args) {
    	System.out.println(gammaQ(20.0,200.0));
    }
}
