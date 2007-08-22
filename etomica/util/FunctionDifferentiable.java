package etomica.util;

public interface FunctionDifferentiable extends Function {

    public double df(int n, double x);

}