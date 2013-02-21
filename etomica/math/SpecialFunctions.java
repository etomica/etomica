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
    public static long factorial(int n){
        if(n < 0){
            throw new IllegalArgumentException("Argument less than zero: "+n);
        }
        if (n < 2) {
            return 1;
        }
        long product = 2;
        for (int i=3; i<n+1; i++) {
            product *= i;
        }
        return product;
    }
    
    public static double lnFactorial(int n) {
        if(n < 0) {
            throw new IllegalArgumentException("Argument less than zero: "+n);
        }
        if (n < 2) {
            return 0;
        }
        long p = 1;
        double sum = 0;
        for (int i=n; i>1; i--) {
            p *= i;
            if (p > Integer.MAX_VALUE/(p-1)) {
                sum += Math.log(p);
                p = 1;
            }
        }
        sum += Math.log(p);
        return sum;
    }

    /**
     * The sign function, returning -1 if x < 0, zero if x == 0, and +1 if x > 0.
     */
    public static double sgn(double x) {
        return (x < 0.0) ? -1.0 : ((x > 0.0) ? +1.0 : 0.0);
    }

    /**
     * Returns the ln(gamma), the natural logarithm of the gamma function.
     * Lanczos approximation, with precision ~15 digits
     * coefficients from GSL (GNU Scientific Library) with g=7
     */
    public static double lnGamma(double x) {
        if (x < 0) {
            throw new IllegalArgumentException("x must not be negative");
        }
        if (x == 0) {
            return 0;
        }
        if (x < 0.5) {
            return Math.log(Math.PI / (Math.sin(Math.PI*x))) - lnGamma(1-x);
        }
        double tmp = x + 7 - 0.5;
        double ser = 0.99999999999980993;
        for(int i=7; i>-1; i--) {
            ser += lnGammaCoeff[i]/(x+i);
        }
        return lnsqrt2Pi + (x-0.5)*Math.log(tmp) - tmp + Math.log(ser);
    }
    private static final double lnsqrt2Pi = 0.5*Math.log(2*Math.PI);
    private static final double[] lnGammaCoeff = new double[] {
                                 676.5203681218851, -1259.1392167224028,
            771.32342877765313, -176.61502916214059, 12.507343278686905,
            -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};

    public static double gamma(double x) {
        if (x == 0) {
            return 1;
        }
        if (x < 0.5) {
            return Math.PI / (Math.sin(Math.PI*x)*gamma(1-x));
        }
        return Math.exp(lnGamma(x));
    }
    
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
    	}
    	return gcf(a, x);
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
        double c = 1.0 / Double.MIN_VALUE;
        double d = 1.0 / b;
        double h = d;
        for (int i = 0; i <= imax; i++) {
            double an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (Math.abs(d) < Double.MIN_VALUE)
                d = Double.MIN_VALUE;
            c = b + an / c;
            if (Math.abs(c) < Double.MIN_VALUE)
                c = Double.MIN_VALUE;
            d = 1.0 / d;
            double del = d * c;
            h *= del;
            if (Math.abs(del - 1.0) < epsilon)
                break;
        }
        return Math.exp(-x + a * Math.log(x) - lnGamma(a)) * h;
    }

    /**
     * The confluent hypergeometric function, usually given the symbol 1F1.  This is a porting to Java
     * of the function hyperg_1F1_luke in the file specfunc/hyperg_1F1.c from the GNU Scientific Library.
     * For algorithm, see [Luke, Algorithms for the Computation of Mathematical Functions, p.182]
     * (GSL license info to be added) GSL home page is http://www.gnu.org/software/gsl/
     */
    public static double confluentHypergeometric1F1(double a, double c, double xin) {
        final double RECUR_BIG = 1.0e+50;
        final double GSL_DBL_EPSILON = 2.2204460492503131e-16;
        final int nmax = 5000;
        int n = 3;
        final double x = -xin;
        final double x3 = x * x * x;
        final double t0 = a / c;
        final double t1 = (a + 1.0) / (2.0 * c);
        final double t2 = (a + 2.0) / (2.0 * (c + 1.0));
        double F = 1.0;
        double prec;

        double Bnm3 = 1.0; /* B0 */
        double Bnm2 = 1.0 + t1 * x; /* B1 */
        double Bnm1 = 1.0 + t2 * x * (1.0 + t1 / 3.0 * x); /* B2 */

        double Anm3 = 1.0; /* A0 */
        double Anm2 = Bnm2 - t0 * x; /* A1 */
        double Anm1 = Bnm1 - t0 * (1.0 + t2 * x) * x + t0 * t1
                * (c / (c + 1.0)) * x * x; /* A2 */

        while (true) {
            double npam1 = n + a - 1;
            double npcm1 = n + c - 1;
            double npam2 = n + a - 2;
            double npcm2 = n + c - 2;
            double tnm1 = 2 * n - 1;
            double tnm3 = 2 * n - 3;
            double tnm5 = 2 * n - 5;
            double F1 = (n - a - 2) / (2 * tnm3 * npcm1);
            double F2 = (n + a) * npam1 / (4 * tnm1 * tnm3 * npcm2 * npcm1);
            double F3 = -npam2 * npam1 * (n - a - 2)
                    / (8 * tnm3 * tnm3 * tnm5 * (n + c - 3) * npcm2 * npcm1);
            double E = -npam1 * (n - c - 1) / (2 * tnm3 * npcm2 * npcm1);

            double An = (1.0 + F1 * x) * Anm1 + (E + F2 * x) * x * Anm2 + F3
                    * x3 * Anm3;
            double Bn = (1.0 + F1 * x) * Bnm1 + (E + F2 * x) * x * Bnm2 + F3
                    * x3 * Bnm3;
            double r = An / Bn;

            prec = Math.abs((F - r) / F);
            F = r;

            if (prec < GSL_DBL_EPSILON || n > nmax)
                break;

            if (Math.abs(An) > RECUR_BIG || Math.abs(Bn) > RECUR_BIG) {
                An /= RECUR_BIG;
                Bn /= RECUR_BIG;
                Anm1 /= RECUR_BIG;
                Bnm1 /= RECUR_BIG;
                Anm2 /= RECUR_BIG;
                Bnm2 /= RECUR_BIG;
                Anm3 /= RECUR_BIG;
                Bnm3 /= RECUR_BIG;
            } else if (Math.abs(An) < 1.0 / RECUR_BIG
                    || Math.abs(Bn) < 1.0 / RECUR_BIG) {
                An *= RECUR_BIG;
                Bn *= RECUR_BIG;
                Anm1 *= RECUR_BIG;
                Bnm1 *= RECUR_BIG;
                Anm2 *= RECUR_BIG;
                Bnm2 *= RECUR_BIG;
                Anm3 *= RECUR_BIG;
                Bnm3 *= RECUR_BIG;
            }

            n++;
            Bnm3 = Bnm2;
            Bnm2 = Bnm1;
            Bnm1 = Bn;
            Anm3 = Anm2;
            Anm2 = Anm1;
            Anm1 = An;
        }

        return F;
        //        result->err  = 2.0 * fabs(F * prec);
        //        result->err += 2.0 * GSL_DBL_EPSILON * (n-1.0) * fabs(F);
    }

    protected static double[] lp = new double[]{1};
    protected static double lastLPX = 0;
    protected static int maxLPN = 0;

    /**
     * Calculates Pn(x), Legendre polynomial of order n for the given input value.
     * The method method makes use of a recursive method to calculate the return
     * value, so calling the method for all n in a sequence should not be much
     * more expensive than the call for the highest order.
     */
    public static double calcLegendrePolynomial(int n, double x) {
        if (n==0 || (x == lastLPX && n<=maxLPN)) {
            return lp[n];
        }
        if (x != lastLPX) {
            // x changed.  pretend we have a cached value, but only for n=0
            maxLPN = 0;
        }
        if (n >= lp.length) {
            // increase size by 2 or 50%, whichever is greater
            int newSize = n+2;
            if (newSize < lp.length*3/2) {
                newSize = lp.length*3/2;
            }
            double[] newLP = new double[newSize];
            System.arraycopy(lp, 0, newLP, 0, maxLPN+1);
            lp = newLP;
        }
        lp[1] = x;
        if (maxLPN == 0) maxLPN=1;
        for (int i=maxLPN+1; i<=n; i++) {
            lp[i] = ((2*i-1)*x*lp[i-1] - (i-1)*lp[i-2])/i;
            System.out.println((2*i-1)*x*lp[i-1]/i + " "+(1-i)*lp[i-2]/i+" "+lp[i]);
        }
        maxLPN = n;
        lastLPX = x;
        return lp[n];
    }

    public static void main(String[] args) {
    	System.out.println(confluentHypergeometric1F1(-0.25,0.5,1.0));
        System.out.println();

        for (int i=-10; i<2; i++) {
            double g = SpecialFunctions.gamma(i-0.5);
            System.out.println((i-0.5)+" "+g); 
        }
        System.out.println(0+" "+SpecialFunctions.gamma(0));
        System.out.println(0.25+" "+SpecialFunctions.gamma(0.25)+" "+Math.exp(SpecialFunctions.lnGamma(0.25)));
        System.out.println(0.5+" "+SpecialFunctions.gamma(0.5)+" "+Math.exp(SpecialFunctions.lnGamma(0.5)));
        System.out.println(1+" "+SpecialFunctions.gamma(1));
        System.out.println();
    	for (int i=2; i<1000; i++) {
    	    double lnfac = SpecialFunctions.lnFactorial(i);
    	    System.out.println(i+" "+lnfac+" "+(SpecialFunctions.lnGamma(i+1)-lnfac)/lnfac);
    	}
    }
}
