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

    private double qiqj;
    private final double alpha;

    private static final double SQRT_PI = Math.sqrt(Math.PI);

    public P2Ewald1Real(double qiqj, double alpha) {
        super(Space3D.getInstance());
        this.qiqj = qiqj;
        this.alpha = alpha;
    }

    public void setQiQj(double newQQ) {
        qiqj = newQQ;
    }

    public double getQiQj() {
        return qiqj;
    }

    @Override
    public double u(double r2) {
        double r = Math.sqrt(r2);
        return qiqj*erfc(alpha*r)/r;
    }

    @Override
    public double du(double r2) {
        double r = Math.sqrt(r2);
        return -qiqj * (2 * SQRT_PI * Math.exp(-alpha * alpha * r2) * alpha + erfc(alpha * r) / r);
    }

    @Override
    public void u012add(double r2, double[] u012) {
        double r = Math.sqrt(r2);
        double u = qiqj * erfc(alpha * r) / r;
        double derfc = -qiqj * 2 / SQRT_PI * Math.exp(-alpha * alpha * r2) * alpha;
        u012[0] += u;
        u012[1] += derfc - u;
        u012[2] += -derfc * (2 + alpha * alpha * 2 * r2) + 2 * u;
    }

    @Override
    public double d2u(double r2) {
        double r = Math.sqrt(r2);
        double u = qiqj * erfc(alpha * r) / r;
        double derfc = -qiqj * 2 / SQRT_PI * Math.exp(-alpha * alpha * r2) * alpha;
        return -derfc * (2 + alpha * alpha * 2 * r2) + 2 * u;
    }

    @Override
    public double integral(double rC) {
        return 0;
    }
}
