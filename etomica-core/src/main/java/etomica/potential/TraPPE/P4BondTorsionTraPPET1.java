package etomica.potential.TraPPE;

import etomica.potential.IPotentialBondTorsion;

public class P4BondTorsionTraPPET1 implements IPotentialBondTorsion {

    protected double c0, c1, c2, c3;
    public P4BondTorsionTraPPET1(double c0, double c1, double c2, double c3) {
        this.c0 = c0;
        this.c1 = c1;
        this.c2 = c2;
        this.c3 = c3;
        // System.out.println(vphi +" "+ n +" "+ phi0);
    }

    @Override
    //virial
    public double u(double costheta) {
        double cosSQtheta = costheta * costheta;
        double cosCBtheta = cosSQtheta * costheta;
        return c0 + c1 * (1 + costheta) + c2 * (2 - 2 * cosSQtheta) + c3 * (1 + 4 * cosCBtheta - 3 * costheta);
    }

    public void udu(double costheta, double[] u, double[] du) {
        // cos(2phi) = 2*cosphi*cosphi - 1
        double cosSQtheta = costheta * costheta;
        double cosCBtheta = cosSQtheta * costheta;

        // a0 + a1*(1+cos(phi)) + a2*(1-cos2phi) + a3*(1+cos3phi)
        u[0] = c0 + c1 * (1 + costheta) + c2 * (2 - 2 * cosSQtheta) + c3 * (1 + 4 * cosCBtheta - 3 * costheta);
        du[0] = 12 * c3 * cosSQtheta - 4 * c2 * costheta + c1 - 3 * c3;
    }
}

