package etomica.mappedDensity;

import etomica.math.function.FunctionDifferentiable;
//uniform external field

public class FunctionUniform implements FunctionDifferentiable {

    private final double L;

    public FunctionUniform(double L) {
        this.L = L;
    }

    @Override
    public double df(int nder, double x) {
        switch (nder) {
            case 0:
                return x + L / 2;
            case 1:
                return 1;
            default:
                return 0;
        }
    }

    @Override
    public double f(double x) {
        return df(0, x);
    }
}
