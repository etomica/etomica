package etomica.potential.UFF;

import etomica.potential.IPotentialBondAngle;
import etomica.potential.IPotentialBondInversion;
import etomica.potential.IPotentialBondTorsion;

public class P4BondInversionUFF implements IPotentialBondInversion {
    protected double c0, c1, c2, kijkl, gamma;
    public P4BondInversionUFF(double c0, double c1, double c2, double kijkl ) {
        this.c0 = c0;
        this.c1 = c1;
        this.c2 = c2;
        this.kijkl = kijkl;

    }

    @Override
    public double u(double cosgamma) {
        double cos2gamma = 2*cosgamma*cosgamma - 1;
        return  kijkl * (c0 + c1 * cosgamma + c2 * cos2gamma);
    }

    @Override
    public void udu(double cosgamma, double[] u, double[] du) {
        double cos2gamma = 2*cosgamma*cosgamma - 1;
        u[0] =  kijkl * (c0 + c1 * cosgamma + c2 * cos2gamma);
        du[0] = 4 * c2 * kijkl * cosgamma - c1 *kijkl;
    }
}
