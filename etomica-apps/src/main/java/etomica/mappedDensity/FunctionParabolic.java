package etomica.mappedDensity;

import etomica.math.SpecialFunctions;
import etomica.math.function.FunctionDifferentiable;
import org.apache.commons.math3.special.Erf;

public class FunctionParabolic implements FunctionDifferentiable {
    double arg = 1;  //Vo
    private final double L;
    private final double temperature;

    public FunctionParabolic(double L, double temperature) {
        this.temperature=temperature;
        this.L = L;
    }

    @Override
    public double df(int nder, double x) {
        switch (nder) {

            //c
            case 0:
                double erfx = 1- SpecialFunctions.erfc(x*Math.sqrt(arg/temperature));
                double erfLby2 = 1- SpecialFunctions.erfc(L*Math.sqrt(arg/temperature)/2);
                //    double erfLm2x = 1- SpecialFunctions.erfc(((L-2*x)/2)*Math.sqrt(arg/temperature));
             //   return (Math.sqrt(temperature*3.14159265359/arg)*((erfL)-(erfLm2x)))/2;
                return  ((erfx+erfLby2)*Math.sqrt(temperature*3.14159265359/arg)/2);

            // p
             case 1:
                return Math.exp(-arg*(x)*(x)/temperature);

            // dp/dz
            case 2:
                return (-2*arg*(x)*Math.exp(-arg*(x)*(x)/temperature)/temperature);
            default:
                throw new RuntimeException("can't do that");
        }
    }

    @Override
    public double f(double x) {
        return df(0, x);
    }

}
