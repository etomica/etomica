package etomica.spin.heisenberg;

import etomica.math.SpecialFunctions;
import etomica.util.numerical.BesselFunction;

public class FunctionCosIntegral {

    /**
     * thetaDot for mean-field mapping of XY model, calculation of correlation function via 2nd free-energy derivative
     * @return thetaDot, with components for derivative with respect to field in x and y directions, respectively
     */
    public static double[] getValue(double x, double b, double theta0) {
        double aC0 = cosInt(x-theta0,b)*BesselFunction.I(1,b)/BesselFunction.I(0,b);
        double C1 = coscosInt(x-theta0,b);
        double S1 = sincosInt(x-theta0,b);
        double cosTheta0 = Math.cos(theta0);
        double sinTheta0 = Math.sin(theta0);
        double A = Math.exp(-b*Math.cos(x-theta0));
        return new double[] {A*(aC0 * cosTheta0 - C1 * cosTheta0 + S1 * sinTheta0),
                A*(aC0 * sinTheta0 - C1 * sinTheta0 + S1 * cosTheta0)};
    }

    /**
     * Integral of Exp(b cos(t)) for t from 0 to x, accurate to within 0.01% for all b > 0 and -2Pi < x < 2Pi.
     */
    public static double cosInt(double x, double b) {
        if(x > 2*Math.PI)
            throw new IllegalArgumentException("Angle must be between -2Pi and + 2Pi");

        // put between 0 and 2Pi
        if(x < 0) return -cosInt(-x, b);

        if(b < 0)
            throw new IllegalArgumentException(("b must be non-negative"));

        if(b <= 1.0) return cosIntSmallb(x, b);

        // put between 0 and Pi
        if(x > Math.PI) return 2*Math.PI*BesselFunction.I(0,b) - cosInt(2*Math.PI - x, b);

        if(x < 0.5*Math.PI) return cosIntLargebA(x, b);

        // if we get here, then Pi/2 < x < Pi
        return cosIntLargebB(x, b);
    }

    /**
     * Integral of sin(t) Exp(b cos(t)) for t from 0 to x, accurate to within 0.01% for all b > 0 and -2Pi < x < 2Pi.
     */
    public static double sincosInt(double x, double b) {
        return (Math.exp(b) - Math.exp(b*Math.cos(x)))/b;
    }


    /**
     * Integral of sin(t) cos(t) Exp(b cos(t)) for t from 0 to x, accurate to within 0.01% for all b > 0 and -2Pi < x < 2Pi.
     */
    public static double sincoscosInt(double x, double b) {
        double cx = Math.cos(x);
        return ((b-1)*Math.exp(b) + (1-b*cx)*Math.exp(b*cx))/(b*b);
    }


    /**
         * Integral of cos(t) Exp(b cos(t)) for t from 0 to x, accurate to within 0.01% for all b > 0 and -2Pi < x < 2Pi.
         */
    public static double coscosInt(double x, double b) {
        if(x > 2*Math.PI)
            throw new IllegalArgumentException("Angle must be between -2Pi and + 2Pi");

        // put between 0 and 2Pi
        if(x < 0) return -coscosInt(-x, b);

        if(b < 0)
            throw new IllegalArgumentException(("b must be non-negative"));

        if(b <= 1.0) return coscosIntSmallb(x, b);

        // put between 0 and Pi
        if(x > Math.PI) return 2*Math.PI*BesselFunction.I(1,b) - coscosInt(2*Math.PI - x, b);

        if(x < 0.5*Math.PI) return coscosIntLargebA(x, b);

        // if we get here, then Pi/2 < x < Pi
        return coscosIntLargebB(x, b);
    }

    /**
     * Integral of cos(t)^2 Exp(b cos(t)) for t from 0 to x.  Evaluates via integration by parts using cosInt and coscosInt.
     * Accurate to within 0.015% for b > 0 and x between 0 and 2pi.
     */
    public static double cos2cosInt(double x, double b) {
        return cosInt(x,b) - (coscosInt(x,b) - Math.sin(x)*Math.exp(b*Math.cos(x))) / b;
    }


    /**
     * Integral of Exp(b cos(t)) from 0 to x, accurate to within 0.01% for b < 1.
     * @param x maximum angle, should be in range -2Pi, 2Pi
     * @param b parameter, should be less than 1
     *
     */
    private static double cosIntSmallb(double x, double b) {
        double sx = Math.sin(x);
        double cx = Math.cos(x);
        double sx2 = sx*sx;
        double cx2 = cx*cx;
        double s2x = 2 * cx * sx;
        double s3x = sx*(3*cx2 - sx2);
        double s4x = 4*cx*sx*(cx2 - sx2);
        double s5x = sx*(5*cx2*cx2 - 10*cx2*sx2 + sx2*sx2);
        double s6x = 2*cx*sx*(3*cx2*cx2 - 10*cx2*sx2 + 3*sx2*sx2);

        return x + b*(sx + b*(0.25*(x + sx*cx) + b*((9*sx + s3x)/72. + b*((12*x + 8*s2x + s4x)/768.
                + b*((150*sx + 25*s3x + 3*s5x)/28800. + b*((60*x + 45*s2x + 9*s4x + s6x)/138240.))))));

    }

    /**
     * Integral of cos(t) Exp(b cos(t)) from 0 to x, accurate to within 0.06% for b < 1.
     * @param x maximum angle, should be in range -2Pi, 2Pi
     * @param b parameter, should be less than 1
     *
     */
    //some duplication of calculation could be removed by combining with cosIntSmallb
    private static double coscosIntSmallb(double x, double b) {
        double sx = Math.sin(x);
        double cx = Math.cos(x);
        double sx2 = sx*sx;
        double cx2 = cx*cx;
        double s2x = 2 * cx * sx;
        double s3x = sx*(3*cx2 - sx2);
        double s4x = 4*cx*sx*(cx2 - sx2);
        double s5x = sx*(5*cx2*cx2 - 10*cx2*sx2 + sx2*sx2);
        double s6x = 2*cx*sx*(3*cx2*cx2 - 10*cx2*sx2 + 3*sx2*sx2);

        double b2 = b*b;
        double b3 = b2*b;
        double b4 = b2*b2;
        double b5 = b4*b;

        return sx + 1./24*b2*(s3x + 9*sx) + (
                b4*(25*s3x + 3*s5x + 150*sx))/5760 + 0.5*b*(cx*sx + x) +
                1./192*b3*(8*s2x + s4x + 12*x) + (
                b5*(45*s2x + 9*s4x + s6x + 60*x))/23040;
    }


    /**
     * Integral of Exp(b cos(t)) from 0 to x, accurate to within 0.01% for b > 1 and -Pi/2 < x < Pi/2.
     */
    private static double cosIntLargebA(double x, double b) {
        double b2 = b*b;
        double b4 = b2*b2;
        double x2 = x*x;
        double x4 = x2*x2;
        double x6 = x4*x2;
        return Math.exp(b)/3628800./b4
                * (Math.exp(-b*x*x/2) * x
                  * (945 + 210*b4*b*x6*(-15 + x2) + 315*b*(600 + x2) +
                     63*b2*(-4050 + 1000*x2 + x4) +
                     9*b2*b*(-50400 - 9450*x2 + 1400*x4 + x6) +
                     b4*x2*(-151200 - 17010*x2 + 1800*x4 + x6)
                    )
                    + 945*(-1 - 200*b + 270*b2 + 480*b2*b + 3840*b4)*Math.sqrt(0.5*Math.PI)
                         * (1-SpecialFunctions.erfc(Math.sqrt(0.5*b)*x))/Math.sqrt(b)
        );
    }

    /**
     * Integral of cos(t) Exp(b cos(t)) from 0 to x, accurate to within 0.05% for b > 1 and -Pi/2 < x < Pi/2.
     */
    private static double coscosIntLargebA(double x, double b) {
        double b2 = b*b;
        double b3 = b2*b;
        double b4 = b2*b2;
        double b5 = b4*b;
        double b6 = b3*b3;
        double sqrtb = Math.sqrt(b);
        double x2 = x*x;
        double x4 = x2*x2;
        double x6 = x4*x2;
        double x8 = x4*x4;
        return 1./(14515200.*b5*sqrtb)*Math.exp(b)*(-2*sqrtb*Math.exp(-b*x2/2)* x *
         (8505 + 2835*b*(466 + x2) + 210*b6*x6*(30 - 17*x2 + x4) +
                189*b2*(-8750 + 2330*x2 + 3*x4) +
                9*b3*(-94500 - 61250*x2 + 9786*x4 + 9*x6) +
                b5*x2*(302400 - 117180*x2 - 14310*x4 + 1378*x6 + x8) +
                9*b4*(-302400 - 31500*x2 - 12250*x4 + 1398*x6 + x8)) +
                945*(9 + 1398*b - 1750*b2 - 900*b3 - 2880*b4 + 7680*b5)*
                Math.sqrt(2*Math.PI)*(1.0-SpecialFunctions.erfc(sqrtb*x/Math.sqrt(2))));
        //Math.sqrt(2*Math.PI)*(1.0-org.apache.commons.math3.special.Erf.erfc(sqrtb*x/Math.sqrt(2))));

    }


    /**
     * Integral of Exp(b cos(t)) from 0 to x, accurate to within 0.01% for b > 1 and Pi/2 < x < Pi.
     */
    private static double cosIntLargebB(double x, double b) {

        double pix = Math.PI - x;
        double pix2 = pix*pix;
        double pix4 = pix2*pix2;
        double pix8 = pix4*pix4;
        double pi2 = Math.PI*Math.PI;
        double pi4 = pi2*pi2;
        double pi3 = pi2*Math.PI;
        double pi8 = pi4*pi4;
        double x2 = x*x;
        double x4 = x2*x2;
        double b2 = b*b;
        double b4 = b2*b2;

        return Math.PI* BesselFunction.I(0, b)
                + (-(1./39916800)*
        Math.exp(-b) * pix*(39916800 + 945*b4*b*pix8*pix2 -
                1050*b4*pix8*(-11 + 3*pi2  - 6*Math.PI*x + 3*x2) +
                15*b2*b*pix4*pix2*(7920 + 147*pi4 - 588*pi3*x -
                1540*x2 + 147*x4 +
                14*pi2*(-110 + 63*x2) + Math.PI*(3080*x - 588*x2*x)) -
                15*b2*pix4*(-66528 + 17*pi4*pi2 - 102*pi3*pi2*x +
                7920*x2 - 462*x4 + 17*x2*x4 - 4*pi3*x*(-462 + 85*x2) +
                3*pi4*(-154 + 85*x2) -
                6*Math.PI*x*(2640 - 308*x2 + 17*x4) +
                3*pi2*(2640 - 924*x2 + 85*x4)) +
                b*pix2*(6652800 + pi8 - 8*pi4*pi3*x - 332640*x2 +
                7920*x4 - 110*x4*x2 + x4*x4 +
                2*pi4*pi2*(-55 + 14*x2) + pi3*pi2*(660*x - 56*x2*x) -
                8*pi3*x*(3960 - 275*x2 + 7*x4) +
                10*pi4*(792 - 165*x2 + 7*x4) +
                2*pi2*(-166320 + 23760*x2 - 825*x4 + 14*x4*x2)
                + Math.PI*(665280*x - 31680*x2*x + 660*x4*x - 8*x4*x2*x))));
    }

    /**
     * Integral of cos(t) Exp(b cos(t)) from 0 to x, accurate to within 0.01% for b > 1 and Pi/2 < x < Pi.
     */
    private static double coscosIntLargebB(double x, double b) {

        double pix = Math.PI - x;
        double pix2 = pix*pix;
        double pix3 = pix2*pix;
        double pix4 = pix2*pix2;
        double pix6 = pix4*pix2;
        double pi2 = Math.PI*Math.PI;
        double pi4 = pi2*pi2;
        double x2 = x*x;
        double x4 = x2*x2;
        double b2 = b*b;
        double b3 = b2*b;
        double b4 = b2*b2;
        double eb = Math.exp(-b);

        return -(1./39916800) *
        eb*pix3*(6652800 +
                pix2*(-332640 + 7920*pix2 - 110*pix4 + pix6 +
                4725*b4*pix6 - 4200*b3*pix4*(-11 + 3*pix2) -
                30*b*(-66528 + 7920*pix2 + pix4*(-462 + 17*pix2)) +
                45*b2*pix2*(7920 +
                7*pix2*(-220 + 21*pix2)))) + (1./39916800)*
        eb*pix*(39916800 +
                b*pix2*(6652800 +
                pix2*(-332640 + 7920*pix2 - 110*pix4 + pix6 +
                945*b4*pix6 - 1050*b3*pix4*(-11 + 3*pix2) -
                15*b*(-66528 + 7920*pix2 + pix4*(-462 + 17*pix2)) +
                15*b2*pix2*(7920 +
                7*pix2*(-220 + 21*pix2))))) + Math.PI*BesselFunction.I(1,b);
    }


    public static void main(String[] args) {
        double b = 10;
        double theta0 = 2*Math.PI;

        // Print results to check against Mathematica
        for(double x = 0; x < 2*Math.PI; x+=0.1) {
           // System.out.println("{" + x + ", " + getValue(x, b,theta0)[0]+"},");
            System.out.println("{" + x + ", " +  cos2cosInt(x, b)+"},");
        }
//        System.out.println(cosIntSmallb(-5.7, .52));
//        System.out.println(cosIntLargebA(1.5, 5.52));
//        System.out.println(cosIntLargebB(2.0, 5.5));
        //System.out.println(coscosIntLargebA(1.2, 9.52));
        //System.out.println(coscosIntLargebB(4.0, 5.5));
        System.out.println(coscosIntSmallb(-5.7, 0.52));
    }
}
