package etomica.mappedDensity.mappedDensityfromlatticesite;

import etomica.math.SpecialFunctions;
import etomica.math.function.FunctionDifferentiable;
import org.apache.commons.math3.special.Erf;

public class FunctionUniform implements FunctionDifferentiable {

    private final double msd;


    public FunctionUniform(double msd) {
        this.msd = msd;
    }


    @Override
    public double df(int nder, double r) {
        switch (nder) {
            //c
            case 0:
                 return r*r*r/3;
            // p
            case 1:
                return 1;
            // dp/dr
            case 2:
                return 0;
            default:
                throw new RuntimeException("can't do that");
        }
    }

    @Override
    public double f(double r) {
        return df(0, r);
    }





}
