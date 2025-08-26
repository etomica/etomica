package etomica.potential.COMPASS;

import etomica.potential.IPotential2;
import etomica.potential.IPotentialBondAngle;
import etomica.potential.IPotentialBondAngleHybrid;
import etomica.potential.P3BondAngleHybrid;

public class BondAngleCOMPASS implements IPotentialBondAngleHybrid {
    private double k1= 100.0;// Spring constant gives a measure of the strength of harmonic interaction
    private final boolean bZero;
    private double b0, theta0;

    public BondAngleCOMPASS(double b0, double theta0, double kt){
        bZero = (b0 == 0.0);
        this.theta0 = theta0;
        this.b0 = b0;
        this.k1 = kt;
    }
    public void u012add(double r, double costheta, double[] u012){
        if(bZero){
            u012[0] = u012[1] = u012[2] =1 ;
        }
        double dx = r - b0;
        double dtheta = Math.acos(costheta) - theta0;
        u012[0] = k1 * dx *dtheta;
    }

    public double u(double r, double costheta) {
        if(bZero ) return 1;
        double dx = r - b0;
        double dtheta = Math.acos(costheta) - theta0;
        return k1 * dx *dtheta;
    }

    public void udu(double costheta, double r12, double[] u, double[] du){
        double dx = r12 - b0;
        double dtheta = Math.acos(costheta) - theta0;
        u[0]= k1 * dx *dtheta;
    }
}
