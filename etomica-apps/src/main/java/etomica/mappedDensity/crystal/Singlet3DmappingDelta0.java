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
     * @param rVec coordinate of sphere
     * @param RVec coordinate of reference lattice site
     * @param deltaVec coordinate of measurement site
     * @param sigma width of Gaussian reference distribution, p(r) = exp(-r^2/(2 sigma^2))
     * @return mapping velocities in x, y, z, directions, respectively
     */
    public static Vector xyzDot(Vector rVec, Vector RVec, Vector deltaVec, double sigma) {

        double sigma2 = sigma*sigma;

        Vector rDir = rVec.makeCopy();
        rDir.ME(deltaVec);

        Vector thetaDir = RVec.makeCopy();//this will be used to store theta direction, but not yet
        thetaDir.ME(deltaVec);

        double r2 = rDir.squared()/sigma2;
        double R2 = thetaDir.squared()/sigma2;
        double rp2 = rVec.Mv1Squared(RVec)/sigma2;
        double rp = Math.sqrt(rp2);
        double r = Math.sqrt(r2);
        double R = Math.sqrt(R2);
        double p = Math.exp(-0.5*rp2);
        double term0 = rVec.dot(RVec)/r; // R * cos(theta)
        double Rsint = Math.sqrt(R2 - term0*term0);// R * sin(theta)
        double term1 = (p*rp*sqrt2OverPi - Erf.erf(rp/sqrt2))/(rp2*rp);
        double term2 = Math.exp(-0.5*R2)/(4*Math.PI*sigma2);

        double Ar = term2 * ( 1/r2 + (r2 - R2 + rp2)*term1/(2*r) );
        double At = term2 * Rsint * term1;

        // now compute theta direction
        thetaDir.XE(rDir);
        thetaDir.XE(rDir);
        thetaDir.normalize();
        rDir.TE(Ar/(r*sigma*p));//divide by (r*sigma) to normalize, divide by p to get rDot
        System.out.println(rDir.squared());
        thetaDir.TE(At/p);

        rDir.PE(thetaDir);
        return rDir;
    }

    public static void main(String[] args) {
        Vector rx = new Vector3D(-.5,.3,2);
        Vector Rx = new Vector3D(1.2,0,0);
        Vector deltax = new Vector3D(.7,1.5,-.1);

        System.out.println(Singlet3DmappingDelta0.xyzDot(rx, Rx,deltax, 1));
    }


}
