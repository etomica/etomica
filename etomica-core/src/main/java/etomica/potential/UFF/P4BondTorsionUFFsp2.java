package etomica.potential.UFF;

import etomica.potential.IPotentialBondTorsion;

public class P4BondTorsionUFFsp2 implements IPotentialBondTorsion {
    protected double vphi, n, phi0;
    public P4BondTorsionUFFsp2(double vphi, double n, double phi0) {
        this.vphi = vphi;
        this.n = n;
        this.phi0 = phi0;

        /*double lambda = 0.1332, n;
        double ri, rj;
        double rbo = -lambda * (ri + rj) * Math.log(n);
        vphi = 5 * Math.sqrt(vsp1 * vsp2) * (1 + 4.18 * Math.log(rbo));
*/
    }

    @Override
    public double u(double cosphi) {
        double phi = Math.acos(cosphi);
        return  0.5* vphi * ( 1 - ( Math. cos(n * phi0) * Math.cos(n * phi) ));
    }

    @Override
    public void udu(double costheta, double[] u, double[] du) {
        throw new RuntimeException("Implement me");
    }

}
