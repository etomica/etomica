package etomica.potential.UFF;

import etomica.potential.IPotentialBondTorsion;

public class P4BondTorsionUFF implements IPotentialBondTorsion {
    protected double vphi, n, phi0;
    public P4BondTorsionUFF(double vphi, double n, double phi0) {
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
    public void udu(double cosphi, double[] u, double[] du) {
//n*sin(nphi)/(Math.sq(1-cos(phi)^2))
        double phi = Math.acos(cosphi);
        u[0] = 0.5* vphi * ( 1 - ( Math. cos(n * phi0) * Math.cos(n * phi) ));
        du[0] = (vphi * n * Math.cos(n * phi) * Math.sin(n * phi)) / ( Math.sin(phi));
    }

}
