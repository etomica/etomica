package etomica.potential.TraPPE;

import etomica.potential.IPotentialBondTorsion;

public class P4BondTorsionTraPPET3 implements IPotentialBondTorsion {
    protected double d0, theta0, theta;

    //d0/2 (spi - spi0)^2
    public P4BondTorsionTraPPET3(double d0, double psi0){
        this.d0 = d0;
        this.theta0 = psi0;
    }

    public double u(double costheta) {
        if(costheta > 1){
            theta  = 0;
        } else if (costheta < -1) {
            theta = Math.PI;
        }

        double theta = Math.acos(costheta);
        return (d0/2) * (theta - theta0) * (theta - theta0);
    }

    public void udu(double costheta, double[] u, double[] du){
        if(costheta > 1){
            theta  = 0;
        } else if (costheta < -1) {
            theta = Math.PI;
        } else {
            theta = Math.acos(costheta);
        }


        u[0] = (d0/2) * (theta - theta0) * (theta - theta0);
        du[0] = 0;
    }
}

