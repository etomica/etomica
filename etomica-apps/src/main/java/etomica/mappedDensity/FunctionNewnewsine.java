package etomica.mappedDensity;

import etomica.math.function.FunctionDifferentiable;

public class FunctionNewnewsine implements FunctionDifferentiable {


    private final int n;
    private final double L;

    public FunctionNewnewsine(int n, double L) {
        this.n = n;
        this.L = L;
    }

    @Override
    public double df(int nder, double x) {
        switch (nder) {
            //c
            case 0:
                return (x + L / 2) * 3 - (2.5*(Math.cos(2 * Math.PI * n * x / L) - Math.cos(Math.PI * n)) / (2 * Math.PI * n / L));
            // p
            case 1:
                return 3 + 2.5*Math.sin(2 * Math.PI * n * x / L);
            // dp/dz
            case 2:
                double arg = 2 * Math.PI * n / L;
                return 2.5*Math.cos(arg * x) * arg;
            default:
                throw new RuntimeException("can't do that");
        }
    }

    @Override
    public double f(double x) {
        return df(0, x);
    }


}
