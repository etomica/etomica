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
    private static final double pi2 = Math.PI * Math.PI;
    //private static final double sqrt2OverPi = Math.sqrt(2.0/Math.PI);
    private static final double sqrtPiOver2 = Math.sqrt(0.5*Math.PI);
    private static final double sqrt2 = Math.sqrt(2.0);
    private static final double twoOver2Pi52 = 2.0/(4*Math.PI*Math.PI*sqrtPi*sqrt2);
    private final Map<IntSet, Double> HnkValues = new HashMap<IntSet, Double>();
    private final double[] SnkValues = new double[62];
    private final double[] CnkValues = new double[62];
    private final double[] dSnkValues = new double[62];
    private final double[] dCnkValues = new double[62];
    private final double[] iSnkValues = new double[62];
    private final double[] iCnkValues = new double[62];
    private final double[] auxValues = new double[62];
    private double lastY = Double.NaN;
    private double lastYaux = Double.NaN;
    private int lastN = -1;

    public Singlet3DmappingDelta0() {}

    /**
     * Returns mapping velocities in x, y, z, directions, in a frame in which the density is
     * measured on the z axis. Actual system must be rotated so that input coordinates are in this frame,
     * and results must be rotated back to laboratory frame.
     * @param nmax Number of terms to include in sums.  5 or so should be ok.
     * @param ri radial coordinate of sphere, relative to the origin (measurement site)
     * @param ti polar angle of sphere, measured from z axis
     * @param phii azimuthal angle of sphere, measured from x axis
     * @param R distance of density measurement site from lattice site
     * @param sigma width of Gaussian reference distribution, p(r) = exp(-r^2/(2 sigma^2))
     * @return mapping velocities in x, y, z, directions, respectively
     */
    public double[] xyzDot(int nmax, int kmax, double ri, double ti, double phii, double R, double sigma) {

        double rihat = ri/sigma;
        double Rhat = R/sigma;

        double coshri = Math.cosh(rihat);
        double sinhri = Math.sinh(rihat);
        double sinh2ri = 2. * sinhri * coshri;
        double cosh2ri = coshri * coshri + sinhri * sinhri;

        double costi = Math.cos(ti);
        double sinti = Math.sin(ti);
        double cos2ti = costi * costi - sinti * sinti;
//        double sin2ti = 2. * sinti * costi;


        double p0 = Math.exp(-0.5 * Rhat * Rhat);
        double expRmr2 = Math.exp(-0.5 * (rihat - Rhat)*(rihat - Rhat));
        double expRpr2 = Math.exp(-0.5 * (rihat + Rhat)*(rihat + Rhat));

        double Ar = 0;
        double At = 0;
        double f = 0;
        double g = 0;
        double irf = 0;
        double irg = 0;
        for(int n=1; n<=nmax; n++) {
            double Rkkfact = 1.0; // R^k / k!
            double rsum = 0.0;
            double tsum = 0.0;
            double isum = 0.0;
            for(int k=0; k<=kmax; k++) {
                double hnkRkkfact = Hnk(n,k) * Rkkfact;
                rsum += hnkRkkfact * dSnk(n,k,rihat);
                tsum += hnkRkkfact * Snk(n,k,rihat);
                isum += hnkRkkfact * iSnk(n,k,rihat);
                Rkkfact *= Rhat/(k+1);
            }
            Ar += Math.cos(n*ti) * rsum / n;
            At += Math.sin(n*ti) * tsum;

            if(n % 2 == 0) { //even n
                double enr = Math.exp(-n*rihat);
                f += -n * p0 / pi2 * ( enr/(n*n-1) - p0/sqrtPi/sqrt2 * tsum);
                irf += -n * p0 / pi2 * ((1 - enr)/n/(n*n-1) - p0/sqrtPi/sqrt2 * isum);
            } else { // odd n
                g += -n * p0 / pi2 * ( 0 - p0/sqrtPi/sqrt2 * tsum);
                irg += -n * p0 / pi2 * ( 0 - p0/sqrtPi/sqrt2 * isum);

            }
        }
        Ar *= -p0*p0 * 2*twoOver2Pi52;
        At *= -p0*p0 * 2*twoOver2Pi52;

        double Atc = (f * Math.cos(ti) + g)/rihat;
        double Arc = -(irf * Math.cos(2*ti) + irg * Math.cos(ti))/(ri * ri * Math.sin(ti));

        // add K0'(ri) term (Ar only)
        double k0P = twoOver2Pi52 * p0 / Rhat * (
                expRmr2 * (1 + sqrtPiOver2 * R * ex2erfc((rihat-Rhat)/sqrt2))
                + expRpr2 * (-1 + sqrtPiOver2 * R * ex2erfc((rihat+Rhat)/sqrt2))
        );
//        System.out.println("k0P: "+k0P);
        Ar += k0P;

        // add delta-even sums
        double lnTerm = Math.log((coshri + costi)/(coshri - costi));
        double tanTerm = Math.atan2(4.*sinhri*sinti, cosh2ri + cos2ti  - 2.);
        double SigmaS = 0.25 * (coshri * sinti * lnTerm - sinhri * costi * tanTerm);
        double SigmaC = 0.25 * (2.0 - sinhri * costi * lnTerm - coshri * sinti * tanTerm);

        Ar += -p0/pi2*SigmaC;
        At += p0/pi2*SigmaS;

        //homogeneous-solution addition to fix meet boundary conditions
        double atc = 0.0;


//        Ar /= rihat * rihat * sinti;
//        At /= -rihat * sinti;
        At = -At;

        Ar -= Arc;
        At -= Atc;

        double p = Math.exp(-0.5*(rihat*rihat + Rhat*Rhat - 2*rihat*Rhat*Math.cos(ti)));
        double riDot = Ar/p;
        double thetaiDot = At/p;


        // thetaiDot here is already times ri sin(thetai)
        // riDot here is already times ri^2 sin(theta)
        double xyDot = riDot / (rihat*rihat) + thetaiDot * costi / sinti;
        double zDot =  riDot * costi / (rihat*rihat*sinti) - thetaiDot;
        return new double[] { xyDot * Math.cos(phii), xyDot * Math.sin(phii), zDot};
    }

    /**
     * The function exp(+x^2) erfc(x). Uses asymptotic expansion for x > 5, and explicit function calls for x < 5.
     * Not well tested for negative x (explicit calls are used, so result should be correct)
     * For x > 5, results are accurate within 10^-9 for n = 10, and 10^-7 for n = 5
     * @param nmax number of terms in asymptotic series.
     * @param x argument of function
     * @return exp(+x^2) erfc(x)
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
     * The function ex2erfc with default argument of nmax = 10
     * @param x function argument
     * @return exp(+x^2) erfc(x)
     */
    public static double ex2erfc(double x) {
        return ex2erfc(10,x);
    }

    /**
     * Returns the integral of (cos(t))^k sin(t) cos(nt) for t from 0 to pi
     */
    public double Hnk(int n, int k) {
        if(k > 61) throw new IllegalArgumentException("Maximum k is 61 for calculation of Hnk");

        if( (n+k) % 2 != 0) return 0.0; //Hnk is zero for n+k odd

        //see if stored value is available
        IntSet key = new IntSet(new int[] {n, k});
        Double valueObj = HnkValues.get(key);
        if(valueObj != null) return valueObj;

        //not previously computed;  compute, store, and return value
        double sum = 0;
        for(int m=0; m<=k; m++) {
            sum += binomial(k,m) * (-1./(m-0.5*(k+n+1)) + 1./(m-0.5*(k+n-1)) - 1./(m-0.5*(k-n+1)) + 1./(m-0.5*(k-n-1)));
        }
        sum /= (1l << k+2);//divide by 2^(k+2)

        HnkValues.put(key, sum);
        return sum;
    }

    /**
     * Returns the difference of integrals:
     *    Int[exp(-x^2/2) x^(k+2) cosh(n(y-x)), {x,0,y}]
     *          - cosh(ny) Lim(L->infinity) Int[exp(-x^2/2) x^(k+2) sinh(n(L-x)), {x,0,L}]/sinh(nL)
     * Uses recursion in k and stores previous values until called with a different n or y
     */
    //should replace recursion with a loop
    public double Cnk(int n, int k, double y) {

        if(haveValue(n, k, y, CnkValues)) return CnkValues[k+2];

        CnkValues[k+2] = - n * Snk(n,k-1,y);

        if(k == -1) CnkValues[k+2] += - Math.exp(-0.5*y*y) ;
        else CnkValues[k+2] += -Math.pow(y,k+1) * Math.exp(-0.5*y*y) + (k+1) * Cnk(n,k-2, y);

        return CnkValues[k+2];
    }

    /**
     * Returns derivative of Cnk with respect to y
     */
    public double dCnk(int n, int k, double y) {

        if(haveValue(n, k, y, dCnkValues)) return dCnkValues[k+2];

        dCnkValues[k+2] = - n * dSnk(n,k-1,y);

        if(k == -1) dCnkValues[k+2] += y * Math.exp(-0.5*y*y) ;
        else dCnkValues[k+2] += -Math.pow(y,k) * Math.exp(-0.5*y*y) * (1 + k - y*y) + (k+1) * dCnk(n,k-2, y);

        return dCnkValues[k+2];
    }

    /**
     * Returns the difference of integrals:
     *    Int[exp(-x^2/2) x^(k+2) sinh(n(y-x)), {x,0,y}]
     *          - cosh(ny) Lim(L->infinity) Int[exp(-x^2/2) x^(k+2) cosh(n(L-x)), {x,0,L}]/sinh(nL)
     * Uses recursion in k and stores previous values until called with a different n or y
     */
    public double Snk(int n, int k, double y) {

        if(haveValue(n, k, y, SnkValues)) return SnkValues[k+2];

        SnkValues[k+2] = - n * Cnk(n,k-1,y);

        if(k == -1) SnkValues[k+2] += - Math.exp(-n*y) ;
        else SnkValues[k+2] += (k+1) * Snk(n,k-2, y);

        return SnkValues[k+2];
    }

    /**
     * Returns derivative of Snk with respect to y
     */
    public double dSnk(int n, int k, double y) {

        if(haveValue(n, k, y, dSnkValues)) return dSnkValues[k+2];

        dSnkValues[k+2] = - n * dCnk(n,k-1,y);

        if(k == -1) dSnkValues[k+2] += n * Math.exp(-n*y) ;
        else dSnkValues[k+2] += (k+1) * dSnk(n,k-2, y);

        return dSnkValues[k+2];
    }

    public double iSnk(int n, int k, double y) {

        if(haveValue(n, k, y, iSnkValues)) return iSnkValues[k+2];

        iSnkValues[k+2] = - n * iCnk(n,k-1,y);

        if(k == -1) iSnkValues[k+2] += (-1 + Math.exp(-n*y))/n ;
        else iSnkValues[k+2] += (k+1) * iSnk(n,k-2, y);

        return iSnkValues[k+2];
    }

    public double iCnk(int n, int k, double y) {

        if(haveValue(n, k, y, iCnkValues)) return iCnkValues[k+2];

        iCnkValues[k+2] = - n * iSnk(n,k-1,y);

        if(k == -1) iCnkValues[k+2] += -sqrtPiOver2 * (1 - erfc(y/sqrt2));
        else iCnkValues[k+2] += -iAux(k,y) + (k+1) * iCnk(n,k-2, y);

        return iCnkValues[k+2];
    }

    public double iAux(int k, double y) {

        if(haveAuxValue(k, y)) return auxValues[k];

        auxValues[k] = -Math.pow(y, k) * Math.exp(-0.5*y*y);

        if(k == 1) auxValues[1] += sqrtPiOver2 * (1 - erfc(y/sqrt2));
        else auxValues[k] += k * iAux(k-2, y);

        return auxValues[k];
    }

    private boolean haveValue(int n, int k, double y, double[] values) {
        if(y != lastY || n != lastN) {
            for (int i = 0; i < SnkValues.length; i++) {
                SnkValues[i] = Double.NaN;
                CnkValues[i] = Double.NaN;
                dSnkValues[i] = Double.NaN;
                dCnkValues[i] = Double.NaN;
                iSnkValues[i] = Double.NaN;
                iCnkValues[i] = Double.NaN;
            }
            lastY = y;
            lastN = n;
        }

        //see if stored value is available
        if(!Double.isNaN(values[k+2])) return true;

        //not previously computed; if termination of recursion, compute values here
        if(k == -2) {
            double ey2 = Math.exp(-0.5*y*y);
            double eny = Math.exp(-n*y);
            int nmax = 20;
            double ex20 = ex2erfc(nmax,n/sqrt2);
            double ey20 = ex2erfc(nmax, y/sqrt2);
            double ex2m = ex2erfc(nmax,(n-y)/sqrt2);
            double ex2p = ex2erfc(nmax,(n+y)/sqrt2);
            CnkValues[0] = sqrtPiOver2 * (-eny * ex20 + 0.5 * ey2 * (ex2m - ex2p));
            SnkValues[0] = -0.5 * sqrtPiOver2 * ey2 * (ex2m + ex2p);
            dCnkValues[0] = ey2 + sqrtPiOver2 * n * (eny * ex20 - 0.5 * ey2 * (ex2m + ex2p));
            dSnkValues[0] = 0.5 * sqrtPiOver2 * n * ey2 * (ex2m - ex2p);
            iCnkValues[0] = -0.5*sqrtPiOver2/n * (-2*eny*ex20 + ey2*(ex2m + ex2p));
            iSnkValues[0] = -sqrtPiOver2/n * (1 + 0.5 * ey2 * (-2*ey20 - ex2m + ex2p));
            return true;
        }

        //not termination of recursion; compute in calling routine
        return false;

    }

    private boolean haveAuxValue(int k, double y) {
        if(y != lastYaux) {
            for (int i = 0; i < SnkValues.length; i++) {
                auxValues[i] = Double.NaN;
            }
            lastYaux = y;
        }

        //see if stored value is available
        if(!Double.isNaN(auxValues[k])) return true;

        //not previously computed; if termination of recursion, compute values here
        if(k == 0) {
            double ey2 = Math.exp(-0.5*y*y);
            auxValues[0] = -ey2 + 1;
            return true;
        }

        //not termination of recursion; compute in calling routine
        return false;

    }

    /**
     * Returns n choose k, n!/k!/(n-k)!
     */
    private static double binomial(int n, int k) {
        if(n-k < k) k = n-k;
        double result = 1;
        for(int i=0; i<k; i++) {
            result *= ((double)(n-i))/(k-i);
        }
        return result;
    }

    public static void main(String[] args) {
        Singlet3DmappingDelta0 s3d = new Singlet3DmappingDelta0();

//        int n = 6;
//        int k = 8;
//        System.out.println(n+" "+k+" "+s3d.Hnk(n,k));

        double rx = .75;
        double thetax = 2.*Math.PI/3.;
        double Rx = 0.8;
        int kmx = 5;
        int nmx = 15;

        s3d.xyzDot(nmx, kmx, rx, thetax,0, Rx,1.0);
//        int n = 2;
//        int k = 3;
//        double y = 2.5;
//        System.out.println(s3d.Snk(n, k, y));
//        System.out.println(s3d.Cnk(n, k, y));
//        System.out.println(s3d.dSnk(n, k, y));
//        System.out.println(s3d.dCnk(n, k, y));

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


}
