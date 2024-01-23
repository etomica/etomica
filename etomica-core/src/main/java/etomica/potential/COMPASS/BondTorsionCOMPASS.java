package etomica.potential.COMPASS;

import etomica.potential.IPotential2;
import etomica.potential.IPotentialBondTorsionHybrid;

public class BondTorsionCOMPASS implements IPotentialBondTorsionHybrid {
    private double k1, k2, k3;// Spring constant gives a measure of the strength of harmonic interaction
    private final boolean bZero;
    private double b0, phi0;

    public BondTorsionCOMPASS(double b0, double phi0, double k1t, double k2t, double k3t){
        bZero = (b0 == 0.0);
        this.phi0 = phi0;
        this.b0 = b0;
        this.k1 = k1t;
        this.k2 = k2t;
        this.k3 = k3t;
    }

    public void u012add(double r, double cosphi, double[] u012){
        if(bZero){
            u012[0] = u012[1] = u012[2] =1 ;
        }
        double phi = Math.acos(cosphi);
        double torsion = (k1 * cosphi) + (k2 * Math.cos(2*phi)) + (k3 * Math.cos(3*phi));
        u012[0] = (r-b0) * torsion;
    }

    public double u(double r, double cosphi) {
        if(bZero) return 1;
        double phi = Math.acos(cosphi);
        double torsion = (k1 * cosphi) + (k2 * Math.cos(2*phi)) + (k3 * Math.cos(3*phi));
        return (r-b0) * torsion;
    }
}
