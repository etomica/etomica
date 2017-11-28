/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math;

import etomica.units.Mole;

/**
 * Static-method library of various functions
 */
public final class SpecialFunctions {

    private static final double lnsqrt2Pi = 0.5 * Math.log(2 * Math.PI);
    private static final double[] lnGammaCoeff = new double[]{
            676.5203681218851, -1259.1392167224028,
            771.32342877765313, -176.61502916214059, 12.507343278686905,
            -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
    protected static double[] lp = new double[]{1};
    protected static double lastLPX = 0;
    protected static int maxLPN = 0;

    private SpecialFunctions() {
    }

    /**
     * Complementary error function, computed using the approximant 7.1.26 of Abramowitz & Stegun.
     * Absolute accuracy to within ~10^-7. If more precision is required,
     * consider using {@link org.apache.commons.math3.special.Erf#erfc}, which is more accurate
     * but substantially slower (~10x - 100x) than the method given here.
     *
     * @param x argument to function, may take values from  Double.NEGATIVE_INFINITY to Double.POSITIVE_INFINITY.
     * @return value of complementary error function.
     */
    public static double erfc(double x) {
        if (x < 0) {
            return 2 - erfc(-x);
        }
        double t = 1.0 / (1.0 + 0.3275911 * x);
        return (t * (
                0.254829592 + t * (
                        -0.284496736 + t * (
                                1.421413741 + t * (
                                        -1.453152027 + 1.061405429 * t))))) * Math.exp(-x * x);
    }

    /**
     * The factorial function, n! Will overflow with any input value greater than 20.
     *
     * @param n an integer greater than or equal to 0
     * @return the factorial of n
     * @throws IllegalArgumentException if n < 0
     */
    public static long factorial(int n) {
        if (n < 0) {
            throw new IllegalArgumentException("Argument less than zero: " + n);
        }
        if (n < 2) {
            return 1;
        }
        long product = 2;
        for (int i = 3; i < n + 1; i++) {
            product *= i;
        }
        return product;
    }

    /**
     * The log of the factorial function, ln(n!). Will not overflow even if n! would be large.
     *
     * @param n an integer greater than or equal to 0
     * @return the natural logarithm of the factorial of n.
     * @throws IllegalArgumentException if n < 0
     */
    public static double lnFactorial(int n) {
        if (n < 0) {
            throw new IllegalArgumentException("Argument less than zero: " + n);
        }
        if (n < 2) {
            return 0;
        }
        long p = 1;
        double sum = 0;
        for (int i = n; i > 1; i--) {
            p *= i;
            if (p > Integer.MAX_VALUE / (i - 1)) {
                sum += Math.log(p);
                p = 1;
            }
        }
        sum += Math.log(p);
        return sum;
    }

    /**
     * Returns the ln(gamma), the natural logarithm of the gamma function.
     * Lanczos approximation, with precision ~15 digits
     * coefficients from GSL (GNU Scientific Library) with g=7
     *
     * @param x a positive number.
     * @return the natural logarithm of the gamma function of x.
     * @throws IllegalArgumentException if x is not positive.
     */
    public static double lnGamma(double x) {
        if (x <= 0) {
            throw new IllegalArgumentException("x must be positive");
        }
        if (x < 0.5) {
            return Math.log(Math.PI / (Math.sin(Math.PI * x))) - lnGamma(1 - x);
        }
        double tmp = x + 7 - 0.5;
        double ser = 0.99999999999980993;
        for (int i = 7; i > -1; i--) {
            ser += lnGammaCoeff[i] / (x + i);
        }
        double y = lnsqrt2Pi + (x - 0.5) * Math.log(tmp) - tmp + Math.log(ser);
        return y;
    }

    /**
     * Returns the value of the gamma function of x.
     *
     * @param x a positive number.
     * @return the gamma function of x.
     * @throws IllegalArgumentException if x is not positive.
     */
    public static double gamma(double x) {
        if (x <= 0) {
            throw new IllegalArgumentException("x must be positive");
        }
        if (x < 0.5) {
            return Math.PI / (Math.sin(Math.PI * x) * gamma(1 - x));
        }
        return Math.exp(lnGamma(x));
    }

    /**
     * The confluent hypergeometric function, usually given the symbol 1F1.
     * For information about the function and use of parameters consult the link:
     * http://mathworld.wolfram.com/ConfluentHypergeometricFunctionoftheFirstKind.html
     */
    public static double confluentHypergeometric1F1(double a, double c, double xin) {
        /*
        takes in the numerator parameter of 1F1 (ap), the denominator parameter of 1F1 (cp) and the value of the argument (z)&
    	returns either the value of the rational approximation of 1F1 when it converges to a given tolerance or,if that given
    	tolerance is not met, will ask for different values of ap, cp & z.
    	Also: A and B	will contain the values of the numerator and denominator polynomials, respectively, for all degrees
    	from 0 to n inclusive where n is the maximum degree for which the values of the polynomials are to be calculated
    	*/


        double x = -xin;
        double[] arrayA = new double[101]; //the numerator of the rational approx.
        double[] arrayB = new double[101]; //the denominator of the rational approx.
        double[] arrayR = new double[101]; //the value of the rational approximations.
        double[] arrayD = new double[101]; //difference in subsequent rational approx.
        int n = 100; //number of iterations,if you change this number must also change condition of the last if statement
        int n1 = n + 1;
        double tolerance = 1E-15;//specified tolerance

        double zero = 0.0;
        double one = 1.0;
        double two = 2.0;
        double three = 3.0;


        // to understand the following code refer to rational approximation of 1F1(ap; cp; -z)

        double ct1 = a * x / c;
        double xn3 = zero;
        double xn1 = two;
        double z2 = x / two;
        double ct2 = z2 / (one + c);
        double xn2 = one;

        arrayA[0] = one;
        arrayB[0] = one;
        arrayB[1] = one + (one + a) * z2 / c;
        arrayA[1] = arrayB[1] - ct1;
        arrayB[2] = one + (two + arrayB[1]) * (two + a) / three * ct2;
        arrayA[2] = arrayB[2] - (one + ct2) * ct1;
        ct1 = three;
        double xn0 = three;
        //for i=3,...,n the values of arrayA[1+i] and arrayB[1+i] are calculated using the recurrence relations below*/

        for (int i = 3; i <= n; i++) {
            //calculation of the multipliers for the recursion
            ct2 = z2 / ct1 / (c + xn1);
            double g1 = one + ct2 * (xn2 - a);
            ct2 = ct2 * (a + xn1) / (c + xn2);
            double g2 = ct2 * ((c - xn1) + (a + xn0) / (ct1 + two) * z2);
            double g3 = ct2 * z2 * z2 / ct1 / (ct1 - two) * (a + xn2) / (c + xn3) * (a - xn2);

            //the recurrance relations for arrayA[i+1] and arrayB[i+1] are as follows

            arrayB[i] = g1 * arrayB[i - 1] + g2 * arrayB[i - 2] + g3 * arrayB[i - 3];
            arrayA[i] = g1 * arrayA[i - 1] + g2 * arrayA[i - 2] + g3 * arrayA[i - 3];

            xn3 = xn2;
            xn2 = xn1;
            xn1 = xn0;
            xn0 = xn0 + one;
            ct1 = ct1 + two;
        }//end for

        arrayD[n1 - 1] = zero;
        for (int j1 = 1; j1 <= n + 1; j1++) {
            arrayR[j1 - 1] = (arrayA[j1 - 1]) / (arrayB[j1 - 1]);//rational approximation of 1f1
            if (j1 > 1) {
                arrayD[j1 - 2] = (arrayR[j1 - 1]) - (arrayR[j1 - 2]);
            }
            if (j1 >= 5 && Math.abs(arrayD[j1 - 2] / arrayR[j1 - 1]) <= tolerance) {
                //checking for convergence of the rational approximation to a given tolerance
                //if that tolerance is met then exit the loop and return the value of the approximation
                return arrayR[j1 - 1];
            }//end if
            //if that tolerance is not met within the given numberof iterations then the program will
            //ask you to check the values entered
            if (j1 == n) {
                throw new RuntimeException("please check your the values a, c & xin");
            }
        }//end for

        return arrayR[n];
    }//end main

    /**
     * Calculates Pn(x), Legendre polynomial of order n for the given input value.
     * The method method makes use of recursion to calculate the return
     * value, so calling the method for all n in a sequence should not be much
     * more expensive than the call for the highest order.
     *
     * @param n order of the polynomial; must be non-negative.
     * @param x argument of the polynomial.
     * @throws IllegalArgumentException if n < 0.
     */
    public static double calcLegendrePolynomial(int n, double x) {
        if (n < 0) {
            throw new IllegalArgumentException("n must be non-negative");
        }

        if (n == 0 || (x == lastLPX && n <= maxLPN)) {
            return lp[n];
        }
        if (x != lastLPX) {
            // x changed.  pretend we have a cached value, but only for n=0
            maxLPN = 0;
        }
        if (n >= lp.length) {
            // increase size by 2 or 50%, whichever is greater
            int newSize = n + 2;
            if (newSize < lp.length * 3 / 2) {
                newSize = lp.length * 3 / 2;
            }
            double[] newLP = new double[newSize];
            System.arraycopy(lp, 0, newLP, 0, maxLPN + 1);
            lp = newLP;
        }
        lp[1] = x;
        if (maxLPN == 0) maxLPN = 1;
        for (int i = maxLPN + 1; i <= n; i++) {
            lp[i] = ((2 * i - 1) * x * lp[i - 1] - (i - 1) * lp[i - 2]) / i;
//            System.out.println((2*i-1)*x*lp[i-1]/i + " "+(1-i)*lp[i-2]/i+" "+lp[i]);
        }
        maxLPN = n;
        lastLPX = x;
        return lp[n];
    }


    /**
     * Calculates the modified(or unmodified) Bessel function of the first kind (I or J).
     * Reference: http://mathworld.wolfram.com/ModifiedBesselFunctionoftheFirstKind.html
     * https://en.wikipedia.org/wiki/Bessel_function#Bessel_functions_of_the_first_kind:_J.CE.B1
     *
     * @param modified true return the modified Bessel function of the first kind I and
     *                 false return the unmodified Bessel function of the first kind J.
     * @param order    of the Bessel function and it's a real number.
     * @param x        the argument of the Bessel function
     * @return modified Bessel function of the first kind to machine precision or
     * unmodified Bessel function of the first kind to the precision of 10E-13.
     */
    private static double bessel(boolean modified, double order, double x) {
        double sum = 0;

        if (order < 0) {
            order *= -1;
        }

        double gammaValue = (order == 0) ? 1 : gamma(order);
        double kTerms = 0;
        double factorial = 1;
        double s1 = Math.pow(x / 2.0, order);

        for (int k = 0; true; k++) {
            factorial *= (k == 0) ? 1 : (modified ? k : -k);
            gammaValue *= (order == 0 && k == 0) ? 1 : (order + k);
            kTerms = s1 * Math.pow((x * x) / 4.0, k) / factorial / gammaValue;
            sum += kTerms;
            if (Math.abs(kTerms / sum) < 1E-15) {
                break;
            }
        }
        return sum;
    }

    /**
     * The modified Bessel function of the first kind, I
     * Reference: http://mathworld.wolfram.com/ModifiedBesselFunctionoftheFirstKind.html
     *
     * @param order of the Bessel function
     * @param x     the argument of the Bessel function
     * @return modified Bessel function of the first kind, to machine precision.
     */
    public static double besselI(double order, double x) {
        return bessel(true, order, x);
    }

    /**
     * The unmodified Bessel function of the first kind, J.
     * https://en.wikipedia.org/wiki/Bessel_function#Bessel_functions_of_the_first_kind:_J.CE.B1
     *
     * @param order of the Bessel function
     * @param x     the argument of the Bessel function
     * @return unmodified Bessel function of the first kind, to the precision of 1E-13.
     */
    public static double besselJ(double order, double x) {
        return bessel(false, order, x);
    }

    public static void main(String[] args) {
        System.out.println(SpecialFunctions.bessel(false, 1.0, 0.2));
        System.exit(2);


        System.out.println(confluentHypergeometric1F1(-0.25, 0.5, 1.0));
        System.out.println();

        for (int i = -10; i < 2; i++) {
            double g = SpecialFunctions.gamma(i - 0.5);
            System.out.println((i - 0.5) + " " + g);
        }
        System.out.println(0 + " " + SpecialFunctions.gamma(0));
        System.out.println(0.25 + " " + SpecialFunctions.gamma(0.25) + " " + Math.exp(SpecialFunctions.lnGamma(0.25)));
        System.out.println(0.5 + " " + SpecialFunctions.gamma(0.5) + " " + Math.exp(SpecialFunctions.lnGamma(0.5)));
        System.out.println(1 + " " + SpecialFunctions.gamma(1));
        System.out.println();
//    	for (int i=2; i<1000; i++) {
//    	    double lnfac = SpecialFunctions.lnFactorial(i);
//    	    System.out.println(i+" "+lnfac+" "+(SpecialFunctions.lnGamma(i+1)-lnfac)/lnfac);
//    	}
//        System.out.println(besselI(0.32,0.25)+" "+besselI(-1.59, 0.25));
//        double ap3 = SpecialFunctions.besselI(3, 200);
//        System.out.println("ap3 = " + ap3);
        double a = 1E-24 / Mole.UNIT.fromSim(1);
        System.out.println(a * a * 3007.044183696613 + " " + a * a * 198.22529144211072);
        System.out.println(a * a * (-28532.603135414935 + 11560.282978358246) + " " + a * a * Math.sqrt(198.22529144211072 * 198.22529144211072 + 33.767924918154755 * 33.767924918154755));

    }
}
