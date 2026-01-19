package etomica.potential.OPLS_AA;

import etomica.potential.IPotentialBondTorsion;

public class P4TorsionOPLS implements IPotentialBondTorsion {

    protected double v1, v2, v3, f1, f2, f3;

    public P4TorsionOPLS(double v1, double v2, double v3, double f1, double f2, double f3){
        this.f1 = f1;
        this.f2 = f2;
        this.f3 = f3;
        this.v1 = v1;
        this.v2 = v2;
        this.v3 = v3;
    }

    @Override
    public double u(double cosphi){
        double phi = Math.acos(cosphi);
        double u = (0.5*v1*(1+Math.cos(phi+f1))) + (0.5*v2*(1-Math.cos(2*phi +f2))) + ((v3/3)*(1+Math.cos(3*phi +f3)));
        return u;
    }

    @Override
    public void udu(double cosphi, double[] u, double[] du){
        double phi = Math.acos(cosphi);
        u[0] = v1/2*(1+Math.cos(phi+f1)) + v2/2*(1-Math.cos(2*phi +f2)) + v3/2*(1+Math.cos(3*phi +f3));
        du[0] = (-v1*Math.sin(phi + f1) + v2*Math.sin(2*phi + f2) - Math.sin(3*phi +f3))/(-Math.sin(phi));
    }
}
