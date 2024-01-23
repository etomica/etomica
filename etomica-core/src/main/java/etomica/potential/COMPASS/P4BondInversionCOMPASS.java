package etomica.potential.COMPASS;

import etomica.potential.IPotential2;
import etomica.potential.IPotentialBondTorsion;

public class P4BondInversionCOMPASS implements IPotential2 {
    private double k1 =100, zhi;

    public P4BondInversionCOMPASS(double k1t){
        this.k1 = k1t;
    }

    public void u012add(double zhi, double[] u012){
        u012[0] = k1 * zhi * zhi;
    }

    public double u (double zhi){
        return k1*zhi*zhi;
    }
}
