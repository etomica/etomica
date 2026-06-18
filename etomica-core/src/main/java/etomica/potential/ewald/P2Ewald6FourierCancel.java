package etomica.potential.ewald;

import etomica.potential.IPotential2;
import etomica.potential.TruncationFactory;

import static etomica.math.SpecialFunctions.factorial;

public class P2Ewald6FourierCancel implements IPotential2 {

    public static IPotential2 makeTruncated(double sigmaI, double epsilonI, double sigmaJ, double epsilonJ, double alpha6, TruncationFactory tf) {
        return tf.make(new P2Ewald6FourierCancel(sigmaI, epsilonI, sigmaJ, epsilonJ, alpha6));
    }

    private final double alpha6Sq;
    private final double alpha66;
    private final double Bij, f6;

    public P2Ewald6FourierCancel(double sigmaI, double epsilonI, double sigmaJ, double epsilonJ, double alpha6) {
        this.alpha6Sq = alpha6 * alpha6;
        this.alpha66 = alpha6Sq * alpha6Sq * alpha6Sq;
        double localBij = 0;
        for (int k = 0; k <= 6; k++) {
            long ck = factorial(6) / (factorial(6 - k) * factorial(k));
            double bik = 0.25 * Math.pow(sigmaI, k) * Math.sqrt(ck * epsilonI);
            double bjk = 0.25 * Math.pow(sigmaJ, 6-k) * Math.sqrt(ck * epsilonJ);
            localBij += bik * bjk;
        }
        this.Bij = localBij;
        f6 = 4*Math.pow(0.5*(sigmaI+sigmaJ),6)*Math.sqrt(epsilonI*epsilonJ);
    }

    @Override
    public double u(double r2) {
        double a2 = r2 * alpha6Sq;
        double a4 = a2*a2;
        return alpha66*(-Bij*(1+a2+a4/2)*Math.exp(-a2) + f6)/(a4*a2);
    }

    @Override
    public double du(double r2) {
        double a2 = r2 * alpha6Sq;
        double a4 = a2*a2;
        double a6 = a4*a2;
        double e = Math.exp(-a2);
        // (2a + 2a3) e/a6 - 2(1 + a2 + a4/2) e/a5 - 6(1 + a2 + a4/2) e/a7
        return alpha66 * (Bij *(6 + 6*a2 + 3*a4 + a6)*e - 6*f6)/a6;
    }

    @Override
    public void u012add(double r2, double[] u012) {
        double a2 = r2 * alpha6Sq;
        double a4 = a2 * a2;
        double a6 = a4 * a2;
        double e = Math.exp(-a2);
        u012[0] += alpha66 * (-Bij * (1 + a2 + a4 / 2) * e + f6)/ a6;
        u012[1] += alpha66 * (Bij * (6 + 6 * a2 + 3 * a4 + a6) * e - 6*f6) / a6;
        u012[2] += alpha66 * (-Bij * (42 + 42 * a2 + 21 * a4 + (7 + 2 * a2) * a6) * e + 42*f6) / a6;
    }

    @Override
    public double d2u(double r2) {
        double a2 = r2 * alpha6Sq;
        double a4 = a2*a2;
        double a6 = a4*a2;
        double e = Math.exp(-a2);
        return alpha66 * (-Bij *(42 + 42*a2 + 21*a4 + (7+2*a2)*a6)*e + 42*f6)/a6;
    }

}
