package etomica.mappedDensity.crystal;

import etomica.space.Vector;
import etomica.space3d.Vector3D;
import org.apache.commons.math3.special.Erf;

/**
 * Returns three dimensional mapping velocity
 * Uses mapping formulation with origin located at delta function (rather than lattice site)
 */
public class Singlet3DmappingDelta0 {

    private static final double sqrt2OverPi = Math.sqrt(2.0/Math.PI);
    private static final double sqrt2 = Math.sqrt(2.0);

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

        double r2 = rVec.squared()/sigma2;
        double R2 = RVec.squared()/sigma2;
        double rp2 = rVec.Mv1Squared(RVec)/sigma2;
        double rp = Math.sqrt(rp2);
        double r = Math.sqrt(r2);
        double rDotR = rVec.dot(RVec)/(r*sigma2); // R * cos(theta)
        double p = Math.exp(-0.5*rp2);
        double Rsint = Math.sqrt(R2 - rDotR*rDotR);// R * sin(theta)
        double term1 = (p*rp*sqrt2OverPi - Erf.erf(rp/sqrt2))/(rp2*rp);
        double term2 = Math.exp(-0.5*R2)/(4*Math.PI*sigma2);

        double Ar = term2 * ( 1/r2 + (r2 - R2 + rp2)*term1/(2*r) );
        double At = term2 * Rsint * term1;

        // minus (theta direction)
        Vector thetaDir = RVec.makeCopy();
        thetaDir.PEa1Tv1(-rDotR/r,rDir); //rDotR = r.R/|r|
//        thetaDir.XE(rDir);
//        thetaDir.XE(rDir);
        thetaDir.normalize();
        rDir.TE(Ar/(r*sigma*p));//divide by (r*sigma) to normalize rDir, divide by p to get rDot
        thetaDir.TE(-At/p);//thetaDir is actually -(theta direction) so need to multiply by -1

        rDir.PE(thetaDir);
        return rDir;
    }

    public static void main(String[] args) {
        Vector rx = new Vector3D(sqrt2,0,sqrt2);
        Vector Rx = new Vector3D(0,0,1.5);

        System.out.println(Singlet3DmappingDelta0.xyzDot(rx, Rx, 1.1));

    }


}
