package etomica.modules.vle;

import etomica.math.SpecialFunctions;
import etomica.statmech.LennardJones;

/**
 * Class that provides method to compute the 2nd virial coefficient for the 
 * Lennard-Jones + point-quadrupole model potential.
 * 
 * @author kofke
 *
 */
public class LJQB2 {

    /**
     * Returns B2 for the LJ+quadrupole model.  All quantities in units of sigma and epsilon.
     * @param T
     * @param Q
     * @return
     */
    public static double B2(double T, double Q) {
        final double gamma712 = 1.528709197087111; //Gamma(7/12)
        final double gamma1312 = 0.9582856821728326; //Gamma(13/12)
        double result = LennardJones.B2(T);
        result += -0.7/T/T*Q*Q * Math.pow(T, 1.0/12.0) / 6.0 / Math.pow(2.0, 1.0/6.0) *
            (Math.sqrt(T) * gamma712 * SpecialFunctions.confluentHypergeometric1F1(7./12., 0.5, 1.0/T)
                    + 2.0 * gamma1312 * SpecialFunctions.confluentHypergeometric1F1(13./12., 1.5, 1.0/T));
        return result;
    }
}
