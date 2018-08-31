package etomica.mappedDensity;

import etomica.math.SpecialFunctions;
import etomica.math.function.FunctionDifferentiable;
import org.apache.commons.math3.special.Erf;

public class FunctionPhiparabolicpsumofgaussians implements FunctionDifferentiable {

    private final double L;
    private final double a1;
    private final double b1;
    private final double c1;
    private final double a2;
    private final double b2;
    private final double c2;
    private final double a3;
    private final double b3;
    private final double c3;
    private final double a4;
    private final double b4;
    private final double c4;

    public FunctionPhiparabolicpsumofgaussians(double L, double a1, double b1, double c1,double a2, double b2, double c2,double a3, double b3, double c3,double a4, double b4, double c4) {
         this.L = L;
        this.a1 = a1;
        this.b1 = b1;
        this.c1 = c1;
        this.a2 = a2;
        this.b2 = b2;
        this.c2 = c2;
        this.a3 = a3;
        this.b3 = b3;
        this.c3 = c3;
        this.a4 = a4;
        this.b4 = b4;
        this.c4 = c4;
    }


    @Override
    public double df(int nder, double x) {
        switch (nder) {
            //c
            case 0:
                double b1plusLby2byc1 = 1- SpecialFunctions.erfc( (b1+L/2)/c1);
                double b2plusLby2byc2 = 1- SpecialFunctions.erfc( (b2+L/2)/c2);
                double b3plusLby2byc3 = 1- SpecialFunctions.erfc( (b3+L/2)/c3);
                double b4plusLby2byc4 = 1- SpecialFunctions.erfc( (b4+L/2)/c4);
                double b1minusLby2byc1 = 1- SpecialFunctions.erfc( (b1-L/2)/c1);
                double b2minusLby2byc2 = 1- SpecialFunctions.erfc( (b2-L/2)/c2);
                double b3minusLby2byc3 = 1- SpecialFunctions.erfc( (b3-L/2)/c3);
                double b4minusLby2byc4 = 1- SpecialFunctions.erfc( (b4-L/2)/c4);
                double b1pluszbyc1 = 1- SpecialFunctions.erfc( (b1+x)/c1);
                double b2pluszbyc2 = 1- SpecialFunctions.erfc( (b2+x)/c2);
                double b3pluszbyc3 = 1- SpecialFunctions.erfc( (b3+x)/c3);
                double b4pluszbyc4 = 1- SpecialFunctions.erfc( (b4+x)/c4);
                double b1minuszbyc1 = 1- SpecialFunctions.erfc( (b1-x)/c1);
                double b2minuszbyc2 = 1- SpecialFunctions.erfc( (b2-x)/c2);
                double b3minuszbyc3 = 1- SpecialFunctions.erfc( (b3-x)/c3);
                double b4minuszbyc4 = 1- SpecialFunctions.erfc( (b4-x)/c4);
                return (0.5*Math.sqrt(3.141592653589)*((a1*c1*(b1plusLby2byc1-b1minuszbyc1))+(a2*c2*(b2plusLby2byc2-b2minuszbyc2))+(a3*c3*(b3plusLby2byc3-b3minuszbyc3))+(a4*c4*(b4plusLby2byc4-b4minuszbyc4))+(a1*c1*(-b1minusLby2byc1+b1pluszbyc1))+(a2*c2*(-b2minusLby2byc2+b2pluszbyc2))+(a3*c3*(-b3minusLby2byc3+b3pluszbyc3))+(a4*c4*(-b4minusLby2byc4+b4pluszbyc4))));

            // p
            case 1:
                return (a1*Math.exp(-((x-b1)/c1)*((x-b1)/c1))+a2*Math.exp(-((x-b2)/c2)*((x-b2)/c2))+a3*Math.exp(-((x-b3)/c3)*((x-b3)/c3))+a4*Math.exp(-((x-b4)/c4)*((x-b4)/c4))+a1*Math.exp(-((x+b1)/c1)*((x+b1)/c1))+a2*Math.exp(-((x+b2)/c2)*((x+b2)/c2))+a3*Math.exp(-((x+b3)/c3)*((x+b3)/c3))+a4*Math.exp(-((x+b4)/c4)*((x+b4)/c4)));

            // dp/dz
            case 2:
                return ((-2*(x-b1)*a1*Math.exp(-((x-b1)/c1)*((x-b1)/c1))/(c1*c1))-(2*(x-b2)*a2*Math.exp(-((x-b2)/c2)*((x-b2)/c2))/(c2*c2))-(2*(x-b3)*a3*Math.exp(-((x-b3)/c3)*((x-b3)/c3))/(c3*c3))-(2*(x-b4)*a4*Math.exp(-((x-b4)/c4)*((x-b4)/c4))/(c4*c4))-(2*(x+b1)*a1*Math.exp(-((x+b1)/c1)*((x+b1)/c1))/(c1*c1))-(2*(x+b2)*a2*Math.exp(-((x+b2)/c2)*((x+b2)/c2))/(c2*c2))-(2*(x+b3)*a3*Math.exp(-((x+b3)/c3)*((x+b3)/c3))/(c3*c3))-(2*(x+b4)*a4*Math.exp(-((x+b4)/c4)*((x+b4)/c4))/(c4*c4)));

            default:
                throw new RuntimeException("can't do that");
        }
    }

    @Override
    public double f(double x) {
        return df(0, x);
    }



}
