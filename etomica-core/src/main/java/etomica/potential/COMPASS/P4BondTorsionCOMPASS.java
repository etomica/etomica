package etomica.potential.COMPASS;

import etomica.potential.IPotentialBondTorsion;

public class P4BondTorsionCOMPASS implements IPotentialBondTorsion {

    protected double k1, k2, k3;

    public P4BondTorsionCOMPASS(double k1t, double k2t, double k3t){
        this.k1 = k1t;
        this.k2 = k2t;
        this.k3 = k3t;
    }

    @Override
    public double u(double cosphi) {
        double phi = Math.acos(cosphi);
        return  k1* ( 1-Math.cos(phi)) + k2*(1-Math.cos(2*phi)) + k3 *(1-Math.cos(3*phi));
    }

    @Override
    public void udu(double cosphi, double[] u, double[] du) {
        double phi = Math.acos(cosphi);
        if(cosphi > 1){
            phi  = 0;
        } else {
            phi = Math.PI;
        }
        u[0] = k1* ( 1-Math.cos(phi)) + k2*(1-Math.cos(2*phi)) + k3 *(1-Math.cos(3*phi));
        du[0] = k1 * Math.sin(phi) + k2 * Math.sin(2*phi) + k3* Math.sin(3*phi);
    }
}
