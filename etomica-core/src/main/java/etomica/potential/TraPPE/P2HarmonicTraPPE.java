package etomica.potential.TraPPE;

import etomica.potential.IPotential2;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;

public class P2HarmonicTraPPE implements IPotential2 {
    private double k = 100.0;
    private final boolean r0Zero;
    private double r0;
    public P2HarmonicTraPPE(double k) {
        this(k, 0.0);
    }

    public P2HarmonicTraPPE(double w, double r0) {
        setSpringConstant(w);
        r0Zero = (r0 == 0.0);
        setR0(r0);
    }
    public P2HarmonicTraPPE(boolean r0Zero) {
        this.r0Zero = r0Zero;
    }

    public void u012add(double r2, double[] u012) {
        if (r0Zero) {
            u012[0] = u012[1] = u012[2] = k*r2/2;
            u012[0] *= 1;
            return;
        }
        double r = Math.sqrt(r2);
        // Unit kjoulepmole = new UnitRatio(new PrefixedUnit(Prefix.KILO, Joule.UNIT), Mole.UNIT);
        // double wnew = kjoulepmole.fromSim(w);
        double dx = r - r0;
        //double denergy = wnew*dx*dx;
        // System.out.println( " Energy Bonding : " + denergy  +" " +  dx );

        u012[0] = k*dx*dx;
        // u012[1] = 2*w*r*dx;
        //u012[2] = 2*w*r2;
    }

    public double u(double r2) {
        if(r0Zero) return k*r2;
        double dx = Math.sqrt(r2) - r0;
        return k*dx*dx;
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        if(r0Zero) return k*r2/2;
        double r = Math.sqrt(r2);
        return k*r*(r-r0)/2;
    }

    public double d2u(double r2) {
        return k*r2/2;
    }

    public double getSpringConstant() {return k;}

    public void setSpringConstant(double factor) {
        k = factor;
    }

    public Dimension getSpringConstantDimension() {
        return new CompoundDimension(new Dimension[]{Energy.DIMENSION, Length.DIMENSION},new double[]{1,-2});
    }

    public double getR0() {
        return r0;
    }

    public void setR0(double r0) {
        this.r0 = r0;
    }

    public Dimension getR0Dimension() {
        return Length.DIMENSION;
    }

}
