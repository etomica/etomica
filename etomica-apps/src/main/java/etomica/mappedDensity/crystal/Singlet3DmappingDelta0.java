package etomica.mappedDensity.crystal;

import etomica.math.SpecialFunctions;
import etomica.space.Vector;
import etomica.space3d.Vector3D;

import static org.apache.commons.math3.special.Erf.erfc;

/**
 * Returns three dimensional mapping velocity
 * Uses mapping formulation with origin located at delta function (rather than lattice site)
 */
public class Singlet3DmappingDelta0 {

    private static final double sqrt2OverPi = Math.sqrt(2.0/Math.PI);
    private static final double sqrt2 = Math.sqrt(2.0);
    private static final double sqrtPi = Math.sqrt(Math.PI);
    private static final double twoPi32 = Math.pow(2.*Math.PI,1.5);
    private static final double sqrtPiOver2 = Math.sqrt(0.5*Math.PI);


    public Singlet3DmappingDelta0() {}

    /**
     * Returns mapping velocities in x, y, z, directions
     * @param rVec coordinate of sphere relative to measurement site
     * @param RVec coordinate of reference lattice site relative to measurement site
     * @param sigma width of Gaussian reference distribution, p(r) = exp(-r^2/(2 sigma^2))
     * @return mapping velocities in x, y, z, directions, respectively
     */

    //rVec is x - a
    //RVec is R - a

    public static Vector xyzDot(Vector rVec, Vector RVec, double sigma) {
        double sigma2 = sigma*sigma;

        double s = Math.sqrt(rVec.squared()); // |x - a|

        if (s < 1e-12) {
            throw new IllegalArgumentException("Mapping velocity is singular at rVec = 0");
        }

        Vector n = rVec.makeCopy();
        n.normalize();

        double q = twoPi32 * sigma2*sigma;

        // p(a), where RVec = R - a, so |a - R| = |RVec|
        double R2 = RVec.squared()/sigma2;
        double pa = Math.exp(-0.5*R2);

        // d = a - R = -RVec
        Vector d = RVec.makeCopy();
        d.TE(-1);

        double b = n.dot(d);
        double b2 = b*b;

        double z = (s + b)/(sqrt2*sigma);

        double bracket =
                sigma2 * (s - b)
                        + sigma * sqrtPiOver2 * (sigma2 + b2) * ex2erfc(z);

        n.TE(pa/(q*s*s) * bracket);
        return n;
    }

    /**
     * Returns mapping velocities in x, y, z, directions, excluding delta-function contribution.
     * Must be multiplied by beta * p(delta-R)
     * @param rtVec coordinate of sphere relative to its lattice site
     * @param rt magnitude of rtVec, |rtVec|; distance ot sphere from its lattice site
     * @param sigma width of Gaussian reference distribution, p(r) = exp(-r^2/(2 sigma^2))
     * @return mapping velocities in x, y, z, directions, respectively
     */
    public static Vector xyzDotA(Vector rtVec, double rt, double sigma) {
        Vector rtDir = rtVec.makeCopy();
        rtDir.TE(1/rt);//normalize

        double term = (sqrt2OverPi * rt/sigma + ex2erfc(rt/(sqrt2*sigma)))/(4*Math.PI * rt*rt);
        rtDir.TE(term);
        return rtDir;
    }

    /**
     * Returns mapping velocities in x, y, z, directions, including only delta-function contribution.
     * Must be multiplied by beta * p(delta-R)
     * @param rVec coordinate of sphere relative to its measurement site, ri - r
     * @param r magnitude of rVec, |rVec|; distance ot sphere from its lattice site
     * @return mapping velocities in x, y, z, directions, respectively
     */
    public static Vector xyzDotB(Vector rVec, double r) {
        Vector rDir = rVec.makeCopy();

        double term = 1/(4 * Math.PI * r*r*r);
        rDir.TE(term);
        return rDir;
    }


    public static void main(String[] args) {
        Vector rx = new Vector3D(1.2*sqrt2,0,sqrt2);
        Vector Rx = new Vector3D(0,0,1.5);

        System.out.println(Singlet3DmappingDelta0.xyzDot(rx, Rx, 1.1));
//        (0.03407535602662447, 0.0, 0.01268644143574632)

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




    /* Functions below here are for the density distributions using a Gaussian approximation for the delta function */


    /**
     * Returns the 3D mapping velocity for a Gaussian-smoothed delta source.
     *
     * This replaces delta(x-a) by a normalized Gaussian of width eps:
     *
     *     G_eps(x-a) = exp[-|x-a|^2/(2 eps^2)] / ((2 pi)^(3/2) eps^3)
     *
     * Coordinate convention:
     *   rVec = x - a  x = location of atom, a = measurement site
     *   RVec = R - a  R = lattice site for reference density distribution
     *
     * The returned vector is v_eps(x; a, R), satisfying
     *
     *   div[ p0(x) v_eps(x) ]
     *     = p0(x) G_eps(x-a) - P_eps(a) p0(x)/q
     *
     * where p0(x) = exp[-|x-R|^2/(2 sigma^2)].
     *
     * For eps -> 0 this approaches the point-source mapping, but for finite eps
     * it is nonsingular near x = a.
     *
     * @param rVec coordinate of atom relative to measurement site, x - a
     * @param RVec coordinate of lattice site relative to measurement site, R - a
     * @param sigma width of Gaussian reference p0
     * @param eps width of Gaussian-smoothed delta source
     * @return mapping velocity vector
     */
    public static Vector xyzDotGaussian(Vector rVec, Vector RVec, double sigma, double eps) {
        if (sigma <= 0.0) {
            throw new IllegalArgumentException("sigma must be positive");
        }
        if (eps <= 0.0) {
            throw new IllegalArgumentException("eps must be positive");
        }

        double sigma2 = sigma*sigma;
        double eps2 = eps*eps;
        double sum2 = sigma2 + eps2;

        /*
         * P_eps(a) = integral p0(x) G_eps(x-a) dx
         *
         * Since RVec = R - a, |a - R|^2 = |RVec|^2.
         */
        double P =
                Math.pow(sigma2/sum2, 1.5)
                        * Math.exp(-0.5*RVec.squared()/sum2);

        /*
         * mu = (eps^2 R + sigma^2 a)/(sigma^2 + eps^2).
         * In coordinates relative to a:
         *
         *   mu - a = eps^2 (R - a)/(sigma^2 + eps^2)
         */
        Vector muRel = RVec.makeCopy();
        muRel.TE(eps2/sum2);

        /*
         * tau^2 = sigma^2 eps^2 / (sigma^2 + eps^2)
         */
        double tau = sigma*eps/Math.sqrt(sum2);

        /*
         * y1 = x - mu = (x-a) - (mu-a)
         */
        Vector y1 = rVec.makeCopy();
        y1.ME(muRel);

        /*
         * y2 = x - R = (x-a) - (R-a)
         */
        Vector y2 = rVec.makeCopy();
        y2.ME(RVec);

        /*
         * p0(x) = exp[-|x-R|^2/(2 sigma^2)]
         */
        double p0x = Math.exp(-0.5*y2.squared()/sigma2);

        /*
         * w = p0 v = P [ E_tau(x-mu) - E_sigma(x-R) ]
         */
        Vector w = gaussianFieldE(y1, tau);
        Vector eSigma = gaussianFieldE(y2, sigma);
        w.ME(eSigma);
        w.TE(P/p0x);

        return w;
    }

    /**
     * Electric-field-like vector generated by a normalized 3D Gaussian charge cloud
     * of width alpha centered at the origin:
     *
     *   div E_alpha(y) = G_alpha(y)
     *
     * where
     *
     *   G_alpha(y) = exp[-|y|^2/(2 alpha^2)] / ((2 pi)^(3/2) alpha^3)
     *
     * and
     *
     *   E_alpha(y)
     *     = y/(4 pi |y|^3)
     *       [ erf(|y|/(sqrt(2) alpha))
     *         - sqrt(2/pi) |y|/alpha exp(-|y|^2/(2 alpha^2)) ].
     *
     * The limiting value at y = 0 is used to avoid 0/0.
     */
    private static Vector gaussianFieldE(Vector y, double alpha) {
        double r2 = y.squared();
        Vector e = y.makeCopy();

        if (r2 == 0.0) {
            e.TE(0.0);
            return e;
        }

        double r = Math.sqrt(r2);

        /*
         * For small r,
         *
         *   E_alpha(y) ~ G_alpha(0) y / 3
         *
         * where G_alpha(0) = 1 / ((2 pi)^(3/2) alpha^3).
         */
        if (r < 1e-8*alpha) {
            double coeff = 1.0/(3.0 * twoPi32 * alpha*alpha*alpha);
            e.TE(coeff);
            return e;
        }

        double z = r/(sqrt2*alpha);
        double expMinusZ2 = Math.exp(-z*z);

        /*
         * erf(z) = 1 - erfc(z)
         *        = 1 - exp(-z^2) ex2erfc(z)
         */
        double erfz = 1.0 - expMinusZ2 * ex2erfc(z);

        double bracket =
                erfz - sqrt2OverPi * (r/alpha) * expMinusZ2;

        double coeff = bracket/(4.0*Math.PI*r2*r);

        e.TE(coeff);
        return e;
    }

    public static double rho0Gaussian(double R2, double sigma, double eps) {
        double sigma2 = sigma*sigma;
        double eps2 = eps*eps;
        double sum2 = sigma2 + eps2;

        double P =
                Math.pow(sigma2/sum2, 1.5)
                        * Math.exp(-0.5*R2/sum2);

        double q = twoPi32 * sigma2*sigma;

        return P/q;
    }
}


