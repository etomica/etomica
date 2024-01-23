package etomica.potential.COMPASS;

import etomica.potential.IPotential2;
import etomica.potential.P2BondBond;

public class BondBondCOMPASS implements P2BondBond {
    private double k1= 100.0;// Spring constant gives a measure of the strength of harmonic interaction
    private final boolean bZero, b1Zero;
    private double b0, b10;

    public BondBondCOMPASS(double b0, double b10, double kt){
        bZero = (b0 == 0.0);
        b1Zero = (b10 == 0.0);
        setR0(b0, b10);
        setSpringConstant(kt);
    }

    public void u012add(double r, double r1, double[] u012){
        if(bZero || b1Zero){
            u012[0] = u012[1] = u012[2] =1 ;
        }
        double dx = r - b0;
        double dx1 = r1 - b10;
        u012[0] = k1 * dx *dx1;
    }

    public double u(double r, double r1, double[] u012) {
        if(bZero || b1Zero) return 1;
        double dx = r - b0;
        double dx1 = r1 - b10;
        return k1 * dx *dx1;
    }



    public void setSpringConstant(double kt){
        k1 =kt;
    }
    public void setR0(double b0, double b10){
        this.b0 = b0;
        this.b10 = b10;
    }
}
