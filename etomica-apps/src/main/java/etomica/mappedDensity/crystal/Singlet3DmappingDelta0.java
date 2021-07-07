package etomica.mappedDensity.crystal;

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

    public Singlet3DmappingDelta0() {}

    /**
     * Returns mapping velocities in x, y, z, directions
     * @param rVec coordinate of sphere relative to measurement site
     * @param RVec coordinate of reference lattice site relative to measurement site
     * @param sigma width of Gaussian reference distribution, p(r) = exp(-r^2/(2 sigma^2))
     * @return mapping velocities in x, y, z, directions, respectively
     */
    public static Vector xyzDot(Vector rVec, Vector RVec, double sigma) {
        double sigma2 = sigma*sigma;
        Vector rDir = rVec.makeCopy();

        double r2 = rVec.squared()/sigma2; // ((r-delta)/sigma)^2
        double R2 = RVec.squared()/sigma2; // ((R-delta)/sigma)^2
        double rp2 = rVec.Mv1Squared(RVec)/sigma2; // ((r - R)/sigma)^2
        double rp = Math.sqrt(rp2); // |r - R|/sigma
        double r = Math.sqrt(r2);   // |r - delta|/sigma
        double rDotR = rVec.dot(RVec)/(r*sigma2); // R * cos(theta)
        double p = Math.exp(-0.5*rp2);
        double Rsint = Math.sqrt(R2 - rDotR*rDotR);// R * sin(theta)
        //double term1 = (p*rp*sqrt2OverPi - Erf.erf(rp/sqrt2))/(rp2*rp);
        double term1 = (rp*sqrt2OverPi + ex2erfc(rp/sqrt2))/(rp2*rp);//not dividing by p here is needed to get rdot
        double term2 = Math.exp(-0.5*R2)/(4*Math.PI*sigma2);

        double Ar = term2 * ( 1/(/*p*/r2) + (r2 - R2 + rp2)*term1/(2*r) );
        double At = term2 * Rsint * term1;
        //if(r2<0.1) System.out.println(r2);

        // minus (theta direction)
        Vector thetaDir = RVec.makeCopy();
        thetaDir.PEa1Tv1(-rDotR/r,rDir); //rDotR = r.R/|r|
//        thetaDir.XE(rDir);
//        thetaDir.XE(rDir);
        thetaDir.normalize();
        rDir.TE(Ar/(r*sigma));//divide by (r*sigma) to normalize rDir
        thetaDir.TE(-At);//thetaDir is actually -(theta direction) so need to multiply by -1

        rDir.PE(thetaDir);
        return rDir;
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

}
