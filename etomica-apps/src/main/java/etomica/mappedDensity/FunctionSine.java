package etomica.mappedDensity;

import etomica.math.function.FunctionDifferentiable;

public class FunctionSine implements FunctionDifferentiable {

    private final int n;
    private final double L;

    public FunctionSine(int n, double L) {
        this.n = n;
        this.L = L;
    }

    @Override
       public double df(int nder, double x) {
           switch (nder) {
               //c
               case 0:
                   return (x + L / 2) * 2 - (Math.cos(2 * Math.PI * n * x / L) - Math.cos(Math.PI * n)) / (2 * Math.PI * n / L);
               // p
               case 1:
                   return 2 + Math.sin(2 * Math.PI * n * x / L);
               // dp/dz
               case 2:
                   double arg = 2 * Math.PI * n / L;
                   return Math.cos(arg * x) * arg;
               default:
                   throw new RuntimeException("can't do that");
           }
       }

//    @Override
//    public double df(int nder, double x) {
//        switch (nder) {

            //c
//            case 0:
//                return ((0.125*L)+(0.25*x)-(0.00404954*L*Math.cos(5*Math.PI*x/L)*Math.cos(5*Math.PI*x/L))+(0.000191145*L*Math.sin(2*10*Math.PI*x/L)));
            // p
//            case 1:
//                return (0.25 + (0.06361*Math.sin(10*Math.PI*x/L)) + (0.01201*Math.sin(2*10*Math.PI*x/L + Math.PI/2)));
             // dp/dz
//            case 2:
//                double arg = 2 * Math.PI * n / L;
//                return  (((1.99837*Math.cos(10*Math.PI*x/L))-(0.754611*Math.sin(2*10*Math.PI*x/L)))/L);
//            default:
//                throw new RuntimeException("can't do that");
//        }
//    }

    @Override
    public double f(double x) {
        return df(0, x);
    }
}
