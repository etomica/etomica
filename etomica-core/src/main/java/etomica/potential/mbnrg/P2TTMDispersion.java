package etomica.potential.mbnrg;

import etomica.exception.MethodNotImplementedException;
import etomica.potential.IPotential2;

/**
 * 2B pairwise Thole-Type (TTM) interatomic dispersion energy
 * parameters: delta (b), C6 (interatomic dispersion coefficients) u = f(Dr)*C6 / r^6
 */


public class P2TTMDispersion implements IPotential2 {
//    public static IPotential2 makeTruncated(double delta, double c6, TruncationFactory tf) {
//        return tf.make(new P2dispersion(delta, c6));
//    }

    private double delta;
    private double c6;

    public P2TTMDispersion(double delta, double c6){
        setDelta(delta);
        setC6(c6);
    }


    double tang_toennies(double x) {
    /**
     * @brief Calculates the tang toennies damping function
     *
     * Given x=delta*r, calculates the tang-toennies damping function of order n
     * @param[in] Order of the TT function, n = 6
     * @param[in] x Value of the variable. It is defined as x=delta*r, where delta is the b constant from the TTM potential,
     * and r is the distance between the two atoms involved
     * @return The value of the damping function
     */

        double one_over_6 = 1.0 / 6.0;
        double one_over_5 = 0.2;
        double one_over_4 = 0.25;
        double one_over_3 = one_over_6 * 2.0;
        double one_over_2 = 0.5;

        double sum6 = 1.0 + x * one_over_6;
        double sum5 = 1.0 + sum6 * x * one_over_5;
        double sum4 = 1.0 + sum5 * x * one_over_4;
        double sum3 = 1.0 + sum4 * x * one_over_3;
        double sum2 = 1.0 + sum3 * x * one_over_2;
        double sum = 1.0 + sum2 * x;

        return 1.0 - sum * Math.exp(-x);
    }
    /**
     * The energy u
     */

    public double u(double r2){
        double r = Math.sqrt(r2);
        double r6 = r2 * r2 * r2;
        double tt = tang_toennies(delta*r);

        return tt * c6 / r6;
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {

        throw new MethodNotImplementedException();

    }
    /**
     * The second derivative of the pair energy, times the square of the
     * separation: r^2 d^2u/dr^2.
     */
    public double d2u(double r2) {

        throw new MethodNotImplementedException();
    }



        /**
         * Integral used for corrections to potential truncation.
         */
//    public double integral(Space space, double rC) {
//        throw new MethodNotImplementedException("Integral for long-range correction for Exp-6 not yet implemented");
//    }

    private void setC6(double c6) {
        this.c6 = c6;
    }
    public double getDelta() {
        return delta;
    }
    public double getC6() {
        return c6;
    }

    private void setDelta(double delta) {
        this.delta = delta;
    }


}
