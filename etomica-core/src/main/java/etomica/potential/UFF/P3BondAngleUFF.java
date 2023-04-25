package etomica.potential.UFF;

import etomica.potential.IPotentialBondAngle;

public class P3BondAngleUFF implements IPotentialBondAngle {
    protected double c0, c1, c2, kijk;
    public P3BondAngleUFF( double c0, double c1, double c2, double kijk) {
        this.c0 = c0;
        this.c1 = c1;
        this.c2 = c2;
        this.kijk = kijk;
    }

    @Override
    public double u(double costheta) {
        double cos2theta = 2 * costheta * costheta - 1;
        return kijk * ( c0 + c1 * costheta + c2 * cos2theta);
    }

    @Override
    public void udu(double costheta, double[] u, double[] du) {
        double cos2theta = 2 * costheta * costheta - 1;
        u[0] =  kijk * ( c0 + c1 * costheta + c2 * cos2theta);
        du[0] =  kijk*(c1 + 4*costheta*c2);
    }
}
