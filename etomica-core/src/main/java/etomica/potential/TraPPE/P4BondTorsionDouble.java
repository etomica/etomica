package etomica.potential.TraPPE;


import etomica.potential.IPotentialBondTorsion;

public class P4BondTorsionDouble implements IPotentialBondTorsion {
    protected double d0, phi0;
    public P4BondTorsionDouble(double d0, double phi0) {
        this.d0 = d0;
        this.phi0 = phi0;
        // System.out.println(vphi +" "+ n +" "+ phi0);
    }

    @Override
    //virial
    public double u(double cosphi) {
        double phi = Math.acos(cosphi);
        return 0.5*d0*(phi-phi0)*(phi-phi0);
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

        u[0] = 0.5*d0*(phi-phi0)*(phi-phi0);
        //System.out.println("Energy Torsion: "+u[0] );
        //du[0] = (vphi * n * Math.cos(n * phi) * Math.sin(n * phi)) / ( Math.sin(phi));
        du[0] = 0.5*d0*(phi-phi0)*(phi-phi0);
    }
}

