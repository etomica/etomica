package etomica.potential.GAFF;

import etomica.potential.IPotentialBondTorsion;

public class P4BondTorsion2 implements IPotentialBondTorsion {

    protected final double k1, k2;
    protected final int n1, n2;
    protected final double d1, d2;

    public P4BondTorsion2(double k1, int n1, double delta1Degrees,
                        double k2, int n2, double delta2Degrees) {
        this.k1 = k1;
        this.n1 = n1;
        this.d1 = (delta1Degrees);

        this.k2 = k2;
        this.n2 = n2;
        this.d2 = (delta2Degrees);
    }

    @Override
    public double u(double costheta) {
        double theta = Math.acos(clamp(costheta));
        return k1 * (1.0 + Math.cos(n1 * theta - d1))
                + k2 * (1.0 + Math.cos(n2 * theta - d2));
    }

    @Override
    public void udu(double costheta, double[] u, double[] du) {
        double c = clamp(costheta);
        double theta = Math.acos(c);

        double arg1 = n1 * theta - d1;
        double arg2 = n2 * theta - d2;

        u[0] = k1 * (1.0 + Math.cos(arg1))
                + k2 * (1.0 + Math.cos(arg2));

        double dUdTheta = -k1 * n1 * Math.sin(arg1)
                - k2 * n2 * Math.sin(arg2);

        double sinTheta = Math.sqrt(Math.max(0.0, 1.0 - c * c));
        du[0] = sinTheta < 1e-12 ? 0.0 : dUdTheta * (-1.0 / sinTheta);
    }

    protected double clamp(double x) {
        if (x > 1.0) return 1.0;
        if (x < -1.0) return -1.0;
        return x;
    }
}
