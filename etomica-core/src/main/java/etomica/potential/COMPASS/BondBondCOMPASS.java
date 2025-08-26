package etomica.potential.COMPASS;

import etomica.potential.P2BondBond;

public class BondBondCOMPASS implements P2BondBond {
    private double k1= 100.0;// Spring constant gives a measure of the strength of harmonic interaction
    private final boolean r2Zero, r1Zero;
    private double r20, r10;

    public BondBondCOMPASS(double b20, double b10, double kt){
        r1Zero = (b10 == 0.0);
        r2Zero = (b20 == 0.0);
        setR0(b10, b20 );
        setSpringConstant(kt);
    }

    public void u012add( double r1,double r2, double[] u012){
        if( r1Zero || r2Zero){
            u012[0] = u012[1] = u012[2] =1 ;
        }
        double dx1 = r1 - r10;
        double dx2 = r2 - r20;
        u012[0] = k1 * dx2 *dx1;
    }

    public double u(double r1 , double r2 ) {
        double dx2 = r2 - r20;
        double dx1 = r1 - r10;
        return k1 * dx2 *dx1;
    }

    public double du(double r1, double r2 ){
        double dx2 = r2 - r20;
        double dx1 = r1 - r10;
        return 2 * k1 * dx2 *dx1;
    }

    public void setSpringConstant(double kt){
        k1 =kt;
    }
    public void setR0(double b10,double b20 ){
        this.r10 = b10;
        this.r20 = b20;
    }
}
