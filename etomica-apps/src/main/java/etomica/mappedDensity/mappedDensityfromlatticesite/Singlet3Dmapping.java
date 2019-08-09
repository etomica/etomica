package etomica.mappedDensity.mappedDensityfromlatticesite;
import etomica.math.SpecialFunctions;

import static org.apache.commons.math3.special.Erf.erfc;
public class Singlet3Dmapping {

    private static final double sqrtPi = Math.sqrt(Math.PI);
    private static final double sqrt2OverPi = Math.sqrt(2.0/Math.PI);
    private static final double sqrt2 = Math.sqrt(2.0);

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
        System.out.println(ex2erfc(10, -10.01));

//        double ri = 0.05;
//        double r = 1.3;
//        double sigma = 1.8;
//        double thetai = 2 * Math.PI / 6.;
//
//        double[] result = xyzDot(10,ri, thetai,0,r,sigma);
//        System.out.println(result[0]+", "+result[1]+", "+result[2]);
    }

}
