package etomica.utility;
//Java2 imports
//import java.util.LinkedList;

import java.util.LinkedList;

public class OdeSolver implements java.io.Serializable {
    
    /**
     * Performs Runge-Kutta integration for fixed number of constant-sized steps.
     *
     * @param n number of integration steps
     * @param h integration step size (constant)
     * @param xy0 integration starting point, initial x and y values
     * @param rhs right-hand side of the differential equation being integrated
     */
    public static Variables[] rungeKutta4(int n, double h, Variables xy0, Rhs rhs) {
        int ny = xy0.y.length;
        double[] k1 = new double[ny];
        double[] k2 = new double[ny];
        double[] k3 = new double[ny];
        double[] f;
        Variables[] variables = new Variables[n]; //array that is returned
        variables[0] = xy0;
        int iStep = 1;
        while(iStep < n) {
            Variables xy = new Variables(xy0);
            f = rhs.dydx(xy);
            for(int i=0; i<ny; i++) {
                k1[i] = h*f[i];
                xy.y[i] = xy0.y[i] + 0.5*k1[i];
            }
            xy.x = xy0.x + 0.5*h;
            f = rhs.dydx(xy);
            for(int i=0; i<ny; i++) {
                k2[i] = h*f[i];
                xy.y[i] = xy0.y[i] + 0.5*k2[i];
            }
            f = rhs.dydx(xy);
            for(int i=0; i<ny; i++) {
                k3[i] = h*f[i];
                xy.y[i] = xy0.y[i] + k3[i];
            }
            xy.x = xy0.x + h;
            f = rhs.dydx(xy);
            for(int i=0; i<ny; i++) {
                xy.y[i] = xy0.y[i] + (k1[i] + 2.0*(k2[i] + k3[i]) + h*f[i])/6.0;
            }
            variables[iStep] = xy;
            iStep++;
            xy0 = xy;
        }//end while
        return variables;
    }
    
    /**
     * Runge-Kutta integration using an adaptive step-size.
     *
     * @param xy0 integration starting point, initial x and y values
     * @param xf  integration end point, final x value
     * @param eps desired absolute accuracy in dependent variables
     * @param rhs right-hand side of the differential equation being integrated
     */
    public static Variables[] rungeKuttaAdaptive(Variables xy0, double xf, double eps, Rhs rhs) {
        int ny = xy0.y.length;
        double[][] k = new double[6][ny];
        double[] y = new double[ny];
        double[] ys = new double[ny];
        int n = 100; //initial choice for number of steps
        double h = (xf - xy0.x)/(double)n;
        LinkedList values = new LinkedList();
        values.add(xy0);
        boolean last = false;
        while(true) {
            Variables xy = new Variables(xy0);
            for(int m=0; m<ny; m++) {y[m] = ys[m] = xy.y[m];}
            for(int i=0; i<6; i++) {
                xy.x = xy0.x + a[i]*h;
                for(int m=0; m<ny; m++) {
                    xy.y[m] = xy0.y[m];
                    for(int j=0; j<i; j++) {
                        xy.y[m] += b[i][j]*k[j][m];
                    }
                }
                double[] f = rhs.dydx(xy);
                for(int m=0; m<ny; m++) {
                    k[i][m] = h*f[m];
                    y[m] += c[i]*k[i][m];
                    ys[m] += cs[i]*k[i][m];
                }//end of m-loop
            }//end of i-loop
            double delta = 0.0;
            for(int m=0; m<ny; m++) {
                delta = Math.max(delta, Math.abs(y[m]-ys[m]));
            }
            if(delta > eps) { //decrease h, redo step
                double hnew = 0.9*h*Math.pow(eps/delta,0.20);
                h = (hnew < 0.1*h) ? 0.1*h : hnew;
                last = false;
            }
            else { //step is acceptable, increase h for next step
                xy.x = xy0.x + h;
                for(int m=0; m<ny; m++) {xy.y[m] = y[m];}
                double hnew = 0.9*h*Math.pow(eps/delta,0.25);
                h = (hnew > 5.0*h) ? 5.0*h : hnew;
                xy0 = xy;
                values.add(xy);
                if(xy.x+h > xf) {
                    if(last) break;
                    h = xf - xy.x;
                    last = true;
                }
            }
        }//end while
        return (Variables[])values.toArray(new Variables[values.size()]);
    }
    
    //Cash-Karp parameters
    private static final double[] a = new double[] {0.0, 0.2, 0.3, 0.6, 1.0, 0.875};
    private static final double[][] b = new double[][] {
        {}, {0.2}, {0.075, 0.225}, {0.3, -0.9, 1.2}, {-11./54., 2.5, -70./27., 35./27.},
        {1631./55296., 175./512., 575./13824., 44275./110592., 253./4096.}
    };
    private static final double[] c = new double[] {37./378., 0., 250./621., 125./594., 0., 512./1771.};
    private static final double[] cs = new double[] {2825./27648., 0., 18575./48384., 13525./55296., 277./14336., 0.25};
    /**
     * Test of integration algorithm.  Sample system is 
     *    y1' = -y2;
     *    y2' = +y1;
     * for which the solution is y1 = cos(x), y2 = sin(x).
     * Initial value is x = 0, and integrates to x = 2Pi.
     */
    public static void main(String[] args) {
        
        int n = 101;
        double h = 2.0*Math.PI/(double)(n-1);
        Rhs rhs = new Rhs() {
            double[] f = new double[2];
            public double[] dydx(Variables XY) {
                f[0] = -XY.y[1];  // dy0/dx = -y1
                f[1] = +XY.y[0];  // dy1/dx = +y0
                return f;
            }};
 //       Variables[] xy = rungeKutta4(n, h, new Variables(0.0, new double[] {1.0, 0.0}), rhs);
        Variables[] xy = rungeKuttaAdaptive(
                new Variables(0.0, new double[] {1.0, 0.0}), 
                2*Math.PI, 1.e-7, rhs);
        
        for(int i=0; i<xy.length; i++) {
            System.out.println(xy[i].y[0]+" "+xy[i].y[1]);
            System.out.println(Math.cos(xy[i].x)+" "+Math.sin(xy[i].x)+" "+xy[i].x);
            System.out.println();
        }
    }//end of main
    
    public interface Rhs {
        public double[] dydx(Variables xy);
    }
    public static class Variables implements java.io.Serializable {
        public double x;
        public double[] y;
        public Variables() {this(1);}
        public Variables(int ny) {y = new double[ny];}
        public Variables(double x, double[] y) {this.x = x; this.y = y;}
        public Variables(Variables xy) {
            x = xy.x;
            y = new double[xy.y.length];
            for(int i=0; i<y.length; i++) {y[i] = xy.y[i];}
        }
    }
    
    
}
