package etomica.potential.ewald;

import etomica.potential.Potential2Soft;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.TruncationFactory;
import etomica.space3d.Space3D;

import static etomica.math.SpecialFunctions.erfc;

public class P2Ewald1Real extends Potential2SoftSpherical {

    public static Potential2Soft makeTruncated(double qiqj, double alpha, TruncationFactory tf) {
        return tf.make(new P2Ewald1Real(qiqj, alpha));
    }

    private final double qiqj;
    private final double alpha;

    private static final double SQRT_PI = Math.sqrt(Math.PI);

    public P2Ewald1Real(double qiqj, double alpha) {
        super(Space3D.getInstance());
        this.qiqj = qiqj;
        this.alpha = alpha;
    }

    @Override
    public double du(double r2) {
        double r = Math.sqrt(r2);
        return -qiqj * (2 * SQRT_PI * Math.exp(-alpha*alpha*r2) *alpha + erfc(alpha*r)/r);
    }

    @Override
    public void udu(double r2, double[] u, double[] du) {
        double r = Math.sqrt(r2);
        u[0] += qiqj*erfc(alpha*r)/r;
        du[0] += -qiqj* 2 * SQRT_PI * Math.exp(-alpha*alpha*r2) * alpha;
    }

    @Override
    public double d2u(double r2) {
        double r = Math.sqrt(r2);
        return -qiqj * (2 * SQRT_PI * Math.exp(-alpha*alpha*r2) *(alpha*(1 - alpha*alpha*2*r)) + erfc(alpha*r)/r);
    }

    @Override
    public double integral(double rC) {
        return 0;
    }

    @Override
    public double u(double r2) {
        double r = Math.sqrt(r2);
        return qiqj*erfc(alpha*r)/r;
    }
}
