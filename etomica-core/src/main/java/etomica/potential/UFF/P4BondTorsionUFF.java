package etomica.potential.UFF;

import etomica.potential.IPotentialBondTorsion;
import etomica.units.*;

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
        //System.out.println("Energy Torsion: "+u[0] );
        //du[0] = (vphi * n * Math.cos(n * phi) * Math.sin(n * phi)) / ( Math.sin(phi));
        //du[0] = -vphi * Math.cos(n * phi0);
        if(n==3){
            du[0] =0.5*3*vphi*Math.cos(3*phi0)*Math.sin(3*phi)/Math.sin(phi);
        } else if (n==6) {
            du[0] =0.5*6*vphi*Math.cos(6*phi0)*Math.sin(6*phi)/Math.sin(phi);
        } else  {
            du[0] =0.5*2*vphi*Math.cos(2*phi0)*Math.sin(2*phi)/Math.sin(phi);
        }
        if(phi == 0){
            Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT),Mole.UNIT);
            du[0]=5*kcals.toSim(2.119)*(1+4.18*Math.log(1.09));
        }
        //du[0] = -0.5*vphi*n*(1-( Math. cos(n * phi0) * Math.sin(n * phi) )) / Math.sin(phi);
    }

}
