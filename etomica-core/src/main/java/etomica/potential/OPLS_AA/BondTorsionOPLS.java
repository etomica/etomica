package etomica.potential.OPLS_AA;

import etomica.potential.IPotentialBondTorsion;
import etomica.potential.P4BondTorsion;
import etomica.space.Space;
import etomica.units.*;

public class BondTorsionOPLS implements IPotentialBondTorsion {
    protected double a0, a1, a2, a3;
    public BondTorsionOPLS( double a0, double a1, double a2, double a3) {
        this.a0 = a0;
        this.a1 = a1;
        this.a2 = a2;
        this.a3 = a3;
    }

    public double u(double costheta) {
        double cos2theta = 2 * costheta * costheta - 1;
        double cos3theta = costheta * (2 * cos2theta - 1);
        double cos4theta = 2 * cos2theta * cos2theta - 1;

        return 0.5 * a0 * (1 + costheta) + 0.5 * a1 * (1 - cos2theta) + 0.5 * a2 * (1 + cos3theta) + 0.5 * a3 * (1 - cos4theta);
    }

    public void udu(double costheta, double[] u, double[] du) {
        double cos2theta = 2 * costheta * costheta - 1;
        double cos3theta = costheta * (2 * cos2theta - 1);
        double cos4theta = 2 * cos2theta * cos2theta - 1;
        u[0] = 0.5 * a0 * (1 + costheta) + 0.5 * a1 * (1 - cos2theta) + 0.5 * a2 * (1 + cos3theta) + 0.5 * a3 * (1 - cos4theta);
    }
}
