package etomica.mappedDensity.mappedDensityfromlatticesite;

import etomica.math.SpecialFunctions;
import etomica.math.function.FunctionDifferentiable;
import org.apache.commons.math3.special.Erf;

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
                return ((-1*msd*6*r*Math.exp(-3*r*r/(2*msd)))+(msd*er*Math.sqrt(msd*6*Math.PI)))/18;
            // p
            case 1:
                return (Math.exp(-r*r/(2*msd/3)));
            // dp/dr
            case 2:
                return (-3*r*Math.exp(-3*r*r/(2*msd))/msd);
            default:
                throw new RuntimeException("can't do that");
        }
    }

    @Override
    public double f(double r) {
        return df(0, r);
    }




}
