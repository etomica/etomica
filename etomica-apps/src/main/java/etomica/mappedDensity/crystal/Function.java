package etomica.mappedDensity.crystal;

import etomica.math.SpecialFunctions;
import etomica.math.function.FunctionDifferentiable;

public class Function implements FunctionDifferentiable  {

    private final double msd;

    public Function(double msd) {
        this.msd = msd;
    }


    @Override
    public double df(int nder, double r) {
        switch (nder) {
            //c
            case 0:
                double er = 1- SpecialFunctions.erfc( r*Math.sqrt(3/(2*msd)));
                return (msd*((-6*r*Math.exp(-3*r*r/(2*msd)))+(Math.sqrt(6*Math.PI*msd)*er))/(18*Math.pow((2*Math.PI*msd/3),1.5)));
            // p
            case 1:
                return (Math.exp(-3*r*r/(2*msd)))/(Math.pow((2*Math.PI*msd/3),1.5));
            // dp/dr
            case 2:
                return (-3*r*(Math.exp(-3*r*r/(2*msd)))/(msd*Math.pow((2*Math.PI*msd/3),1.5)));
            default:
                throw new RuntimeException("can't do that");
        }
    }

    @Override
    public double f(double r) {
        return df(0, r);
    }




}
