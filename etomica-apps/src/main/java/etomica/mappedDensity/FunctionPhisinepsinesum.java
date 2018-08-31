package etomica.mappedDensity;

import etomica.math.function.FunctionDifferentiable;

public class FunctionPhisinepsinesum implements FunctionDifferentiable {

    private final double L;
    private final double aa;
    private final double bb;
    private final double cc;

    public FunctionPhisinepsinesum(double L, double aa, double bb, double cc) {
        this.L = L;
        this.aa = aa;
        this.bb = bb;
        this.cc = cc;
    }

    @Override
    public double df(int nder, double x) {
        switch (nder) {
             //c
            case 0:
                return (((10*aa*Math.PI*(L+(2*x)))-(4*bb*L*Math.cos(5*Math.PI*x/L)*Math.cos(5*Math.PI*x/L))+(cc*L*Math.sin(2*10*Math.PI*x/L)))/(20*Math.PI));
            // p
            case 1:
                return (aa + (bb*Math.sin(10*Math.PI*x/L)) + (cc*Math.sin(2*10*Math.PI*x/L + Math.PI/2)));
            // dp/dz
            case 2:
                 return  (((10*bb*Math.PI*Math.cos(10*Math.PI*x/L))-(20*cc*Math.PI*Math.sin(2*10*Math.PI*x/L)))/L);
            default:
                throw new RuntimeException("can't do that");
        }
    }

    @Override
    public double f(double x) {
        return df(0, x);
    }


}
