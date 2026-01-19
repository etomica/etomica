package etomica.potential.GAFF;

import etomica.potential.IPotentialBondTorsion;
import etomica.units.*;

public class P4BondTorsionGAFF implements IPotentialBondTorsion {
    protected double vphi, n, phi, gamma;
    public P4BondTorsionGAFF(double vphi, double gamma, double n){
        this.vphi = vphi;
        this.n = n;
        this.gamma = gamma;;
    }

    @Override
    public double u(double cosphi){
        double phi = Math.acos(cosphi);
        return 0.5 * vphi * ( 1 + Math.cos(n * phi - gamma)) ;
    }

    @Override
    public void udu(double cosphi, double[] u, double[] du) {
        //n*sin(nphi)/(Math.sq(1-cos(phi)^2))
        double phi = Math.acos(cosphi);
        /*if(cosphi > 1){
            phi  = 0;
        } else {
            phi = Math.PI;
        }*/
        Unit kjoulepmole = new UnitRatio(new PrefixedUnit(Prefix.KILO, Joule.UNIT), Mole.UNIT);
        double wnew = kjoulepmole.fromSim(vphi);
        double dx = n * phi - gamma;
        double denergy = wnew*dx*dx;
        System.out.println(dx + " " + n * phi * 180 / Math.PI +" " + gamma * 180 / Math.PI+ " " + wnew + " Energy : " + denergy );
        u[0] = vphi*( 1 + Math.cos(n * phi - gamma)) ;
        System.out.println(u[0] + " in udu \n");
        //du[0] = (vphi * n * Math.cos(n * phi) * Math.sin(n * phi)) / ( Math.sin(phi));
        du[0] = -  n * vphi * Math.sin(n*phi - gamma);
    }
}
