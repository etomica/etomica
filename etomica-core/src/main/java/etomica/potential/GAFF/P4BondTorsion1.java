package etomica.potential.GAFF;

import etomica.potential.IPotentialBondTorsion;

public class P4BondTorsion1 implements IPotentialBondTorsion {

    protected final double k1;
    protected final int n1;
    protected final double d1;

    public P4BondTorsion1(double k1, int n1, double delta1Degrees) {
        this.k1 = k1;
        this.n1 = n1;
        this.d1 = (delta1Degrees);
    }

    @Override
    public double u(double costheta) {
        double theta = Math.acos(clamp(costheta));
        return k1 * (1.0 + Math.cos(n1 * theta - d1));
    }

    @Override
    public void udu(double costheta, double[] u, double[] du) {
        double c = clamp(costheta);
        double theta = Math.acos(c);

        double arg = n1 * theta - d1;
        u[0] = k1 * (1.0 + Math.cos(arg));

        double dUdTheta = -k1 * n1 * Math.sin(arg);
        double sinTheta = Math.sqrt(Math.max(0.0, 1.0 - c * c));
        du[0] = sinTheta < 1e-12 ? 0.0 : dUdTheta * (-1.0 / sinTheta);
    }

    protected double clamp(double x) {
        if (x > 1.0) return 1.0;
        if (x < -1.0) return -1.0;
        return x;
    }
}
