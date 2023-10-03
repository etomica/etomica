package etomica.potential.UFF;

import etomica.potential.IPotential2;
import etomica.units.*;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;

public class P2HarmonicUFF implements IPotential2 {
    private double w = 100.0;// Spring constant gives a measure of the strength of harmonic interaction
    private final boolean r0Zero;
    private double r0;

    public P2HarmonicUFF(double w) {
        this(w, 0.0);
    }

    public P2HarmonicUFF(double w, double r0) {
        setSpringConstant(w);
        r0Zero = (r0 == 0.0);
        setR0(r0);
    }

    public P2HarmonicUFF(boolean r0Zero) {
        this.r0Zero = r0Zero;
    }

    public void u012add(double r2, double[] u012) {
        if (r0Zero) {
            u012[0] = u012[1] = u012[2] = 2*w*r2;
            u012[0] *= 1;
            return;
        }
        double r = Math.sqrt(r2);
        Unit kjoulepmole = new UnitRatio(new PrefixedUnit(Prefix.KILO, Joule.UNIT), Mole.UNIT);
        double wnew = kjoulepmole.fromSim(w);
        double dx = r - r0;
        double denergy = wnew*dx*dx;
       // System.out.println( " Energy Bonding : " + denergy  +" " +  dx );

        u012[0] = w*dx*dx;
        u012[1] = 2*w*r*dx;
        u012[2] = 2*w*r2;
    }

    public double u(double r2) {
        if(r0Zero) return w*r2;
        double dx = Math.sqrt(r2) - r0;
        return w*dx*dx;
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        if(r0Zero) return 2*w*r2;
        double r = Math.sqrt(r2);
        return 2*w*r*(r-r0);
    }

    public double d2u(double r2) {
        return 2*w*r2;
    }

    public double getSpringConstant() {return w;}

    public void setSpringConstant(double factor) {
        w = factor;
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
