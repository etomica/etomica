package etomica.potential.TraPPE;

import etomica.potential.IPotentialBondTorsion;

public class P4BondTorsionTraPPET2 implements IPotentialBondTorsion {
    // c0 + c1 cosphi + c2 cos2phi + c3 cos 3phi
    protected double c0, c1, c2, c3, phi0;
    public P4BondTorsionTraPPET2(double c0, double c1, double c2, double c3) {
        this.c0 = c0;
        this.c1 = c1;
        this.c2 = c2;
        this.c3 = c3;
    }

    public double u (double costheta){
        double cosSQtheta = costheta * costheta;
        double cosCBtheta = cosSQtheta * costheta;
        return c0 + c1 * costheta  + c2 * ( 2 * cosSQtheta - 1) + c3 * (4 * cosCBtheta - 3 * costheta);
    }

    public void udu (double costheta, double[] u, double[] du){
        // cos(2phi) = 2*cosphi*cosphi - 1
        double cosSQtheta = costheta * costheta;
        double cosCBtheta = cosSQtheta * costheta;
        u[0] = c0 + c1 * costheta  + c2 * ( 2 * cosSQtheta - 1) + c3 * (4 * cosCBtheta - 3 * costheta);
        du[0] = (c1 - 3 * c3) + (4 * c2 * costheta) + (12 * c3 * cosSQtheta);
    }
}
