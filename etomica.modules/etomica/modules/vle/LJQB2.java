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
        final double gamma1712 = 0.8864821077509258; //Gamma(17/12)
        final double gamma116 = 0.9406558582567716; //Gamma(11/6)
        final double gamma2312 = 0.9675843510790108; //Gamma(23/12)
        final double gamma73 = 1.190639348758999; //Gamma(7/3)
        final double gamma34 = 1.225416702465178; //Gamma(3/4)
        final double gamma54 = 0.9064024770554771; //Gamma(5/4)
        final double two16 = 1.122462048309373; // 2^(1/6)
        final double two56 = 1.781797436280679; // 2^(5/6)
        final double sqrtPI = Math.sqrt(Math.PI);
        final double Q2 = Q*Q;
        final double Q4 = Q2*Q2;
        final double beta = 1.0/T;
        final double sqrtT = Math.sqrt(T);
        double t1 = LennardJones.B2(T);
        double t2 = -7. * Math.PI * Q4 / Math.pow(T, 23.0/12.0) / 60. / two16 *
            (sqrtT * gamma712 * SpecialFunctions.confluentHypergeometric1F1(7./12., 0.5, beta)
                    + 2.0 * gamma1312 * SpecialFunctions.confluentHypergeometric1F1(13./12., 1.5, beta));
        double t3 = 3. * Math.PI * Q2*Q4 / Math.pow(T, 2.5) / 245. *
            (sqrtT * SpecialFunctions.confluentHypergeometric1F1(1.0, 0.5, beta)
                + sqrtPI * SpecialFunctions.confluentHypergeometric1F1(1.5, 1.5, beta));
//        double t3alt = 3. * Math.PI * Q2*Q2*Q2 / Math.pow(T, 2.5) / 245. * 
//            (sqrtT - Math.exp(beta) * sqrtPI * (-2.0 + SpecialFunctions.erfc(beta)));
        double t4 = -71. * Math.PI * Q4*Q4 / Math.pow(T, 37.0/12.0) / 1960. / two56 * 
            (sqrtT * gamma1712 * SpecialFunctions.confluentHypergeometric1F1(17./12., 0.5, beta)
                    + 2.0 * gamma2312 * SpecialFunctions.confluentHypergeometric1F1(23./12., 1.5, beta));
        double t5 =  16. * Math.pow(2.0, 1.0/3.0) * Math.PI * Q4*Q4*Q2 / Math.pow(T, 11.0/3.0) / 4235. * 
            (sqrtT * gamma116 * SpecialFunctions.confluentHypergeometric1F1(11./6., 0.5, beta)
                + 2.0 * gamma73 * SpecialFunctions.confluentHypergeometric1F1(7./3., 1.5, beta));
        double t6 = -18033. * Math.PI*Math.PI * Q4*Q4*Q4 / Math.pow(T, 4.25) / 146578432. *
            (10. * sqrtT / gamma34 * SpecialFunctions.confluentHypergeometric1F1(2.25, 0.5, beta)
                   + 21.0 / gamma54 * SpecialFunctions.confluentHypergeometric1F1(2.75, 1.5, beta) );
//        System.out.println(SpecialFunctions.erfc(1.0));
//        System.out.println("B2LJ+corrections: "+t1+"\t"+t2+"\t"+t3+"\t"+t4+"\t"+t5+"\t"+t6);
        return t1 + t2 + t3 + t4 + t5 + t6;
    }
    
    public static void main(String[] args) {
        double T = 1.0;
        double[] Q2 = {0.0, 0.1, 0.2, 0.4, 1.0};
        
        for(int i=0; i<Q2.length; i++) {
            System.out.println("Q^2, B2 = " + Q2[i] + "\t" + B2(T,Math.sqrt(Q2[i])));
        }
    }
}
