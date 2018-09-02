package etomica.mappedDensity;

import etomica.math.function.FunctionDifferentiable;

public class FunctionLnparabolic implements FunctionDifferentiable {

    private final double L;
    private final double arg=15;

    public FunctionLnparabolic(double L) {
        this.L = L;
    }

    @Override
    public double df(int nder, double x) {
        switch (nder) {
            //c
            case 0:
                return ((L/2)+x+arg*((L*L*L/24)+(x*x*x/3)));
            // p
            case 1:
                return (1 + arg*x*x);
            // dp/dz
            case 2:
                 return (2*x*arg);
            default:
                throw new RuntimeException("can't do that");
        }
    }


    @Override
    public double f(double x) {
        return df(0, x);
    }


}
