package etomica.potential.COMPASS;

import etomica.potential.IPotential2;

public class P2HarmonicCOMPASS implements IPotential2 {
    private double k2 = 100.0, k3=100, k4 =100;// Spring constant gives a measure of the strength of harmonic interaction
    private final boolean r0Zero;
    private double r0;

    public P2HarmonicCOMPASS(double k2, double k3, double k4) {
        this(k2, k3, k4, 0.0);
    }

    public P2HarmonicCOMPASS(double k2t, double k3t, double k4t, double r0) {
        setSpringConstant(k2t, k3t, k4t);
        r0Zero = (r0 == 0.0);
        setR0(r0);
    }

    public P2HarmonicCOMPASS(boolean r0Zero){
        this.r0Zero = r0Zero;
    }

    public void u012add(double r2, double[] u012){
        double r = Math.sqrt(r2);
        double dx = r - r0;
        u012[0] = k2*dx*dx +k3*dx*dx*dx +k4*dx*dx*dx*dx;
        u012[1] = -2*k2*dx -3*k3*dx*dx -4*k4*dx*dx*dx;
        u012[2] = 2*k2 + 6*k3*dx +12*k4*dx*dx;
    }

    public double u(double r2) {
        if(r0Zero) {
            double r = Math.sqrt(r2);
            return k2*r2 + k3*r2*r + k4*r2*r2;
        }
        double dx = Math.sqrt(r2) - r0;
        return k2*dx*dx +k3*dx*dx*dx +k4*dx*dx*dx*dx;
    }

    public double du(double r2) {
        if(r0Zero) {
            double r = Math.sqrt(r2);
            return  -2 *k2*r-3*k3*r2-4*k4*r2*r;
        }
        double dx = Math.sqrt(r2) - r0;
        return -2*k2*dx -3*k3*dx*dx -4*k4*dx*dx*dx;
    }

    public double d2u(double r2) {
        if(r0Zero) {
            double r = Math.sqrt(r2);
            return 2*k2 + 6*k3*r +12*k4*r2;
        }
        double dx = Math.sqrt(r2) - r0;
        return 2*k2 + 6*k3*dx +12*k4*dx*dx;
    }

    public void setSpringConstant(double  k2t, double k3t, double k4t){
        k2 = k2t;
        k3 = k3t;
        k4 = k4t;
    }

    public double getR0() {
        return r0;
    }

    public void setR0(double r0) {
        this.r0 = r0;
    }
}
