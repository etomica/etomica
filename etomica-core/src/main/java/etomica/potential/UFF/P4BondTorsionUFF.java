package etomica.potential.UFF;

import etomica.potential.IPotentialBondTorsion;

public class P4BondTorsionUFF implements IPotentialBondTorsion {
    protected int n;
    protected double vphi, phi0;
    public P4BondTorsionUFF(double vphi, int n, double phi0) {
        this.vphi = vphi;
        this.n = n;
        this.phi0 = phi0;
       // System.out.println(vphi +" "+ n +" "+ phi0);
    }

    @Override
    //virial
    public double u(double cosphi) {
        double phi = Math.acos(cosphi);
        return  0.5* vphi * ( 1 - ( Math. cos(n * phi0) * Math.cos(n * phi) ));
    }

    @Override
    public void udu(double cosphi, double[] u, double[] du) {
        //n*sin(nphi)/(Math.sq(1-cos(phi)^2))
        double phi = Math.acos(cosphi);
        if(cosphi > 1){
            phi  = 0;
        } else if (cosphi < -1) {
            phi = Math.PI;
        }

        u[0] = 0.5* vphi * ( 1 - ( Math. cos(n * phi0) * Math.cos(n * phi) ));
        System.out.println(u[0] + " in udu");
        //du[0] = (vphi * n * Math.cos(n * phi) * Math.sin(n * phi)) / ( Math.sin(phi));
        du[0] = -vphi * Math.cos(n * phi0);
    }

}
