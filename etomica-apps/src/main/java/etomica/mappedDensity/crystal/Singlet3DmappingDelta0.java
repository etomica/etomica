package etomica.mappedDensity.crystal;

import etomica.virial.IntSet;
import org.apache.commons.math3.special.Erf;

import java.util.HashMap;
import java.util.Map;

import static org.apache.commons.math3.special.Erf.erfc;

/**
 * Returns three dimensional mapping velocity
 * Uses mapping formulation with origin located at delta function (rather than lattice site)
 */
public class Singlet3DmappingDelta0 {

    private static final double sqrtPi = Math.sqrt(Math.PI);
    private static final double sqrt2OverPi = Math.sqrt(2.0/Math.PI);
    private static final double sqrtPiOver2 = Math.sqrt(0.5*Math.PI);
    private static final double sqrt2 = Math.sqrt(2.0);
    private static final Map<IntSet, Double> HnkValues = new HashMap<IntSet, Double>();
    private static final double[] SnkValues = new double[62];
    private static final double[] CnkValues = new double[62];
    private static final double[] dSnkValues = new double[62];
    private static final double[] dCnkValues = new double[62];
    private static double lastY = Double.NaN;
    private static int lastN = -1;
    private static final int nmax = 20;

    /**
     * Returns mapping velocities in x, y, z, directions, in a frame in which the density is
     * measured on the z axis. Actual system must be rotated so that input coordinates are in this frame,
     * and results must be rotated back to laboratory frame.
     * @param nmax Number of terms to include in sums.  5 or so should be ok.
     * @param ri radial coordinate of sphere, relative to its lattice site
     * @param ti polar angle of sphere, measured from z axis
     * @param phii azimuthal angle of sphere, measured from x axis
     * @param r distance of density measurement site from lattice site
     * @param sigma width of Gaussian reference distribution, p(r) = exp(-r^2/(2 sigma^2))
     * @return mapping velocities in x, y, z, directions, respectively
     */
    public static double[] xyzDot(int nmax, double ri, double ti, double phii, double r, double sigma) {

        double sigma2 = sigma*sigma;
        double sqrt2pisigma3 = Math.sqrt(2*Math.PI) * sigma * sigma2;
        double sqrt2sigma = Math.sqrt(2.0) * sigma;

        double expr2 = Math.exp(-0.5*r*r/sigma2);
        double expri2 = Math.exp(+0.5*ri*ri/sigma2);
        double costi = Math.cos(ti);

        double cosh2r = Math.cosh(2*r);
        double cosh4r = Math.cosh(4*r);

        double sinh2ri = Math.sinh(2*ri);
        double cosh2ri = Math.cosh(2*ri);
        double sinh4ri = Math.sinh(4*ri);
        double cosh4ri = Math.cosh(4*ri);

        double sinti = Math.sin(ti);
        double cos2ti = Math.cos(2*ti);
        double sin2ti = Math.sin(2*ti);
        double sin4ti = Math.sin(4*ti);
        double cos4ti = Math.cos(4*ti);
        double cosh2ti = Math.cosh(2*ti);

        // thetaiDot * ri * sin(thetai)
        double term1 = sqrt2OverPi*(-sigma*costi + (-Math.PI + 2*ti + Math.PI*(1+sigma2)*costi)/(Math.PI*sigma));

        double term2 = 0.5 * expri2 * (8 * cosh2r * cosh2ri * sin2ti - 4 * sin4ti);
        term2 /= Math.PI * (1 + cos4ti + cosh4r - 4 * cos2ti * cosh2r * cosh2ri + cosh4ri);

        double thetaDotAnalytic =  expr2  * (term1 + term2);

        double thetaDotSum = 0.0;
        for(int n=1; n<=nmax; n++) {
            thetaDotSum += (1+4*n*n*sigma2)/(-1+4*n*n) * Math.sin(2 * n * ti) *
                    (ex2erfc((-ri+2*n*sigma2)/sqrt2sigma) + ex2erfc((+ri+2*n*sigma2)/sqrt2sigma) - sqrt2OverPi/(n*sigma));
        }
        thetaDotSum *= 2 * expr2 / Math.PI;

        double thetaiDot = thetaDotAnalytic + thetaDotSum;//this is times ri * sin(theta)

        // riDot * ri^2 * sin(thetai)
        term1 = -ri * (Math.PI/6. - ti + (-2 + ti*ti)/Math.PI + sinti)/(sigma*sigma2);

        term2 = expri2 * sqrt2OverPi * (-2 * cos2ti * cosh2r * sinh2ri + sinh4ri);
        term2 /= 1 + cos4ti + cosh4r - 4 * cos2ti * cosh2r * cosh2ri + cosh4ri;

        double rDotAnalytic = sqrt2OverPi * expr2 * (term1 + term2);

        double rDotSum = 0.0;
        for(int n=1; n<=nmax; n++) {
            rDotSum += (1+4*n*n*sigma2)/(-1+4*n*n) * Math.cos(2 * n * ti) *
                    (ex2erfc((-ri+2*n*sigma2)/sqrt2sigma) - ex2erfc((+ri+2*n*sigma2)/sqrt2sigma) - ri/(n*n*sqrt2pisigma3));
        }
        rDotSum *= 2 * expr2 / Math.PI;

        double rDot0 = sqrt2OverPi*ri/sigma + ex2erfc(ri/(sqrt2*sigma));
        if(ri < r) rDot0 -= expri2;
        rDot0 *= 2.0 * expr2 / Math.PI;

        double riDot = rDot0 + rDotAnalytic + rDotSum;


        // thetaiDot here is already times ri sin(thetai)
        // riDot here is already times ri^2 sin(theta)
        double xyDot = riDot / (ri*ri) + thetaiDot * costi / sinti;
        double zDot =  riDot * costi / (ri*ri*sinti) - thetaiDot;
        return new double[] { xyDot * Math.cos(phii), xyDot * Math.sin(phii), zDot};
    }

    /**
     * The function exp(+x^2) erfc(x). Uses asymptotic expansion for x > 5, and explicit function calls for x < 5.
     * Not well tested for negative x (explicit calls are used, so result should be correct)
     * For x > 5, results are accurate within 10^-9 for n = 10, and 10^-7 for n = 5
     * @param nmax number of terms in asymptotic series.
     * @param x argument of function
     * @return
     */
    public static double ex2erfc(int nmax, double x) {

        //if (x < 5.) return Math.exp(x * x) * SpecialFunctions.erfc(x);

        if (x < 5.) return Math.exp(x * x) * erfc(x);//apache version of erfc, slower than SpecialFunctions, but more accurate

        double sum = 1.0;
        double prod = 1.0;
        double x2 = 2 * x * x;

        for (int n = 1; n <= nmax; n++) {
            prod *= -(2 * n - 1) / x2;
            sum += prod;
        }

        return sum / x / sqrtPi;

    }

    /**
     * The function ex2erfc with default argument of n = 10
     * @param x
     * @return
     */
    public static double ex2erfc(double x) {
        return ex2erfc(10,x);
    }

    public static void main(String[] args) {

        System.out.println(deltaSnk(3,13, 0.5));
        System.out.println(deltaCnk(3,13, 0.5));

        //System.out.println(binomial(10,4));
//        System.out.println(ex2erfc(10, -10.01));

//        double ri = 0.05;
//        double r = 1.3;
//        double sigma = 1.8;
//        double thetai = 2 * Math.PI / 6.;
//
//        double[] result = xyzDot(10,ri, thetai,0,r,sigma);
//        System.out.println(result[0]+", "+result[1]+", "+result[2]);
    }

    /**
     * Returns the integral of (cos(t))^k sin(t) cos(nt) for t from 0 to pi
     */
    public static double Hnk(int n, int k) {

        if(k > 61) throw new IllegalArgumentException("Maximum k is 61 for calculation of Hnk");

        if( (n+k) % 2 != 0) return 0.0; //Hnk is zero for n+k odd

        //see if stored value is available
        IntSet key = new IntSet(new int[] {n, k});
        Double valueObj = HnkValues.get(key);
        if(valueObj != null) return valueObj.doubleValue();

        //not previously computed;  compute, store, and return value
        double sum = 0;
        for(int m=0; m<=k; m++) {
            sum += binomial(k,m) * (-1./(m-0.5*(k+n+1)) + 1./(m-0.5*(k+n-1)) - 1./(m-0.5*(k-n+1)) + 1./(m-0.5*(k-n-1)));
        }
        sum /= (1l << k+2);//divide by 2^(k+2)

        HnkValues.put(key, new Double(sum));
        return sum;
    }

//    public static double Snk(int n, int k, double y) {
//
//        if(y != lastY) {
//            SnkValues.clear();
//            CnkValues.clear();
//        }
//
//        //see if stored value is available
//        IntSet key = new IntSet(new int[] {n, k});
//        Double valueObj = SnkValues.get(key);
//        if(valueObj != null) return valueObj.doubleValue();
//
//        //not previously computed;  compute, store, and return value
//        double value = -n * Cnk(n,k-1,y) + (k+1) * Snk(n,k-2,y);
//
//        SnkValues.put(key, new Double(value));
//        return value;
//    }

//    /**
//     * Returns the integral of exp(-x^2/2) x^k sinh(n(y-x)) for x from 0 to y.
//     * Uses recursion in k and stores previous values until called with a different n or y.
//     * Could be improved to better handle large n and k
//     */
//
//    public static double Snk(int n, int k, double y) {
//
//        if(y != lastY || n != lastN) {
//            for (int i = 0; i < SnkValues.length; i++) {
//                SnkValues[i] = Double.NaN;
//                CnkValues[i] = Double.NaN;
//            }
//            lastY = y;
//            lastN = n;
//        }
//
//        //see if stored value is available
//        if(!Double.isNaN(SnkValues[k])) return SnkValues[k];
//
//        //not previously computed;  compute, store, and return value
//        if(k == 0) {
//            double exp1 = Math.exp(0.5*n*(n-2*y));
//            double exp2 = Math.exp(2.0*n*y);
//            double erf1 = Erf.erf(n/sqrt2);
//            double erf2 = Erf.erf((y-n)/sqrt2);
//            double erf3 = Erf.erf((y+n)/sqrt2);
//            double erf13 = Erf.erf((y+n)/sqrt2, n/sqrt2); //erf1-erf3
//            //SnkValues[0] = -0.5 * exp1 * sqrtPiOver2 * ((1+exp2)*erf1 + erf2 - exp2*erf3);
//            SnkValues[0] = -0.5 * exp1 * sqrtPiOver2 * (erf1 + erf2 + exp2*erf13);
//            CnkValues[0] = 0.5 * exp1 * sqrtPiOver2 * (erf1 + erf2 - exp2*erf13);
//            return SnkValues[0];
//        }
//
//        SnkValues[k] = -n * Cnk(n,k-1,y);
//
//        if(k == 1) SnkValues[k] += Math.sinh(n*y);
//        else SnkValues[k] += (k-1) * Snk(n,k-2, y);
//
//        return SnkValues[k];
//    }
//
//    /**
//     * Returns the integral of exp(-x^2/2) x^k cosh(n(y-x)) for x from 0 to y.
//     * Uses recursion in k and stores previous values until called with a different n or y
//     */
//    public static double Cnk(int n, int k, double y) {
//
//        if(y != lastY || n != lastN) {
//            for (int i = 0; i < SnkValues.length; i++) {
//                SnkValues[i] = Double.NaN;
//                CnkValues[i] = Double.NaN;
//            }
//            lastY = y;
//            lastN = n;
//        }
//
//        //see if stored value is available
//        if(!Double.isNaN(CnkValues[k])) return CnkValues[k];
//
//        //not previously computed;  compute, store, and return value
//        if(k == 0) {
//            double exp1 = Math.exp(0.5*n*(n-2*y));
//            double exp2 = Math.exp(2.0*n*y);
//            double erf1 = Erf.erf(n/sqrt2);
//            double erf2 = Erf.erf((y-n)/sqrt2);
//            double erf3 = Erf.erf((y+n)/sqrt2);
//            double erf13 = Erf.erf((y+n)/sqrt2, n/sqrt2); //erf1-erf3
//            SnkValues[0] = -0.5 * exp1 * sqrtPiOver2 * (erf1 + erf2 + exp2*erf13);
//            CnkValues[0] = 0.5 * exp1 * sqrtPiOver2 * (erf1 + erf2 - exp2*erf13);
////            SnkValues[0] = -0.5 * exp1 * sqrtPiOver2 * ((1+exp2)*erf1 + erf2 - exp2*erf3);
////            CnkValues[0] = 0.5 * exp1 * sqrtPiOver2 * (-(-1+exp2)*erf1 + erf2 + exp2*erf3);
//            return CnkValues[0];
//        }
//
//        CnkValues[k] = -Math.pow(y,k-1) * Math.exp(-0.5*y*y) - n * Snk(n,k-1,y);
//
//        if(k == 1) CnkValues[k] += Math.cosh(n * y);
//        else CnkValues[k] += (k-1) * Cnk(n,k-2, y);
//
//        return CnkValues[k];
//    }

    /**
     * Returns the difference of integrals:
     *    Int[exp(-x^2/2) x^(k+2) cosh(n(y-x)), {x,0,y}]
     *          - cosh(ny) Lim(L->infinity) Int[exp(-x^2/2) x^(k+2) sinh(n(L-x)), {x,0,L}]/sinh(nL)
     * Uses recursion in k and stores previous values until called with a different n or y
     */
    //should replace recursion with a loop
    public static double deltaCnk(int n, int k, double y) {
        if(y != lastY || n != lastN) {
            for (int i = 0; i < SnkValues.length; i++) {
                dSnkValues[i] = Double.NaN;
                dCnkValues[i] = Double.NaN;
            }
            lastY = y;
            lastN = n;
        }

        //see if stored value is available
        if(!Double.isNaN(dCnkValues[k+2])) return dCnkValues[k+2];

        //not previously computed;  compute, store, and return value
        if(k == -2) {
            dCnkValues[0] = sqrtPiOver2 * (-Math.exp(-n*y) * ex2erfc(nmax,n/sqrt2)
                    + 0.5 * Math.exp(-0.5*y*y) * (ex2erfc(nmax,(n-y)/sqrt2) - ex2erfc(nmax,(n+y)/sqrt2)));
            dSnkValues[0] = -0.5 * sqrtPiOver2 * (Math.exp(-0.5*y*y) * (ex2erfc(nmax,(n-y)/sqrt2) + ex2erfc(nmax,(n+y)/sqrt2)));
            return dCnkValues[0];
        }

        CnkValues[k+2] = - n * deltaSnk(n,k-1,y);

        if(k == -1) CnkValues[k+2] += - Math.exp(-0.5*y*y) ;
        else CnkValues[k+2] += -Math.pow(y,k+1) * Math.exp(-0.5*y*y) + (k+1) * deltaCnk(n,k-2, y);

        return CnkValues[k+2];
    }

    /**
     * Returns the difference of integrals:
     *    Int[exp(-x^2/2) x^(k+2) sinh(n(y-x)), {x,0,y}]
     *          - cosh(ny) Lim(L->infinity) Int[exp(-x^2/2) x^(k+2) cosh(n(L-x)), {x,0,L}]/sinh(nL)
     * Uses recursion in k and stores previous values until called with a different n or y
     */
    public static double deltaSnk(int n, int k, double y) {

        if(y != lastY || n != lastN) {
            for (int i = 0; i < SnkValues.length; i++) {
                dSnkValues[i] = Double.NaN;
                dCnkValues[i] = Double.NaN;
            }
            lastY = y;
            lastN = n;
        }

        //see if stored value is available
        if(!Double.isNaN(dSnkValues[k+2])) return dSnkValues[k+2];

        //not previously computed;  compute, store, and return value
        if(k == -2) {
            dCnkValues[0] = sqrtPiOver2 * (-Math.exp(-n*y) * ex2erfc(nmax,n/sqrt2)
                    + 0.5 * Math.exp(-0.5*y*y) * (ex2erfc(nmax,(n-y)/sqrt2) - ex2erfc(nmax,(n+y)/sqrt2)));
            dSnkValues[0] = -0.5 * sqrtPiOver2 * (Math.exp(-0.5*y*y) * (ex2erfc(nmax,(n-y)/sqrt2) + ex2erfc(nmax,(n+y)/sqrt2)));
            return dSnkValues[0];
        }

        SnkValues[k+2] = - n * deltaCnk(n,k-1,y);

        if(k == -1) SnkValues[k+2] += - Math.exp(-n*y) ;
        else SnkValues[k+2] += (k+1) * deltaSnk(n,k-2, y);

        return SnkValues[k+2];
    }


    private static double binomial(int n, int k) {
        if(n-k < k) k = n-k;
        double result = 1;
        for(int i=0; i<k; i++) {
            result *= ((double)(n-i))/(k-i);
        }
        return result;
    }

}
