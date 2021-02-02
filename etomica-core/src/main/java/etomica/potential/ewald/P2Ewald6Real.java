package etomica.potential.ewald;

import etomica.potential.Potential2Soft;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.TruncationFactory;
import etomica.space.Space;
import etomica.space3d.Space3D;

import static etomica.math.SpecialFunctions.factorial;

public class P2Ewald6Real extends Potential2SoftSpherical {

    public static Potential2Soft makeTruncated(Space space, double sigmaI, double epsilonI, double sigmaJ, double epsilonJ, double alpha6, TruncationFactory tf) {
        return tf.make(new P2Ewald6Real(sigmaI, epsilonI, sigmaJ, epsilonJ, alpha6));
    }

    private final double alpha6Sq;
    private final double alpha66;
    private final double Bij;

    public P2Ewald6Real(double sigmaI, double epsilonI, double sigmaJ, double epsilonJ, double alpha6) {
        super(Space3D.getInstance());
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
    }

    @Override
    public double du(double r2) {
        double a2 = r2 * alpha6Sq;
        double a4 = a2*a2;
        double a6 = a4*a2;
        double e = Math.exp(-a2);
        // (2a + 2a3) e/a6 - 2(1 + a2 + a4/2) e/a5 - 6(1 + a2 + a4/2) e/a7
        return Bij * alpha66 *(6 + 6*a2 + 3*a4 + a6)*e/a6;
    }

    @Override
    public void u012add(double r2, double[] u012) {
        double a2 = r2 * alpha6Sq;
        double a4 = a2 * a2;
        double a6 = a4 * a2;
        double e = Math.exp(-a2);
        u012[0] += -Bij * alpha66 * (1 + a2 + a4 / 2) * e / a6;
        u012[1] += Bij * alpha66 * (6 + 6 * a2 + 3 * a4 + a6) * e / a6;
        u012[2] += -Bij * alpha66 * (42 + 42 * a2 + 21 * a4 + (7 + 2 * a2) * a6) * e / a6;
    }

    @Override
    public void u01TruncationCorrection(double[] uCorrection, double[] duCorrection) {

    }

    @Override
    public double d2u(double r2) {
        double a2 = r2 * alpha6Sq;
        double a4 = a2*a2;
        double a6 = a4*a2;
        double e = Math.exp(-a2);
        return -Bij * alpha66 *(42 + 42*a2 + 21*a4 + (7+2*a2)*a6)*e/a6;
    }

    @Override
    public double integral(double rC) {
        return 0;
    }

    @Override
    public double u(double r2) {
        double a2 = r2 * alpha6Sq;
        double a4 = a2*a2;
        return -Bij*alpha66*(1+a2+a4/2)*Math.exp(-a2)/(a4*a2);
    }
}
