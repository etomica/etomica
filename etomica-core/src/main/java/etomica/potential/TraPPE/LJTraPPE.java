package etomica.potential.TraPPE;

import etomica.potential.IPotential2;
import etomica.potential.P2LennardJones;
import etomica.potential.TruncationFactory;
import etomica.potential.UFF.LJUFF;
import etomica.space.Space;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;

public class LJTraPPE implements IPotential2 {
    public static IPotential2 makeTruncated(double sigma, double epsilon, TruncationFactory tf) {
        return tf.make(new LJTraPPE(sigma, epsilon));
    }

    public LJTraPPE() {
        this(1.0, 1.0);
    }

    public LJTraPPE(double sigma, double epsilon) {
        super();
        setSigma(sigma);
        setEpsilon(epsilon);
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
        double s2 = sigmaSquared/r2;
        double s6 = s2*s2*s2;
        return epsilon4*s6*(s6 - 1.0);
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        double s2 = sigmaSquared/r2;
        double s6 = s2*s2*s2;
        return -epsilon48*s6*(s6 - 0.5);
    }

    public void u012add(double r2, double[] u012) {
        double s2 = sigmaSquared / r2;
        double s6 = s2 * s2 * s2;
        u012[0] += epsilon4 * s6 * (s6 - 1.0);
        u012[1] += -epsilon48 * s6 * (s6 - 0.5);
        u012[2] += epsilon624 * s6 * (s6 - _168div624);
    }

    /**
     * The second derivative of the pair energy, times the square of the
     * separation:  r^2 d^2u/dr^2.
     */
    public double d2u(double r2) {
        double s2 = sigmaSquared/r2;
        double s6 = s2*s2*s2;
        return epsilon624*s6*(s6 - _168div624);
    }

    /**
     * Integral used for corrections to potential truncation.
     */
    public double integral(Space space, double rC) {
        double A = space.sphereArea(1.0);  //multiplier for differential surface element
        int D = space.D();                 //spatial dimension
        double rc = sigma / rC;
        double rCD = space.powerD(rC);
        double rc3 = rc * rc * rc;
        double rc6 = rc3 * rc3;
        double rc12 = rc6 * rc6;
        return 4.0 * epsilon * rCD * A * (rc12 / (12. - D) - rc6 / (6. - D));
    }

    /**
     * Accessor method for Lennard-Jones size parameter.
     */
    public double getSigma() {return sigma;}
    /**
     * Mutator method for Lennard-Jones size parameter.
     * Does not adjust potential cutoff for change in sigma.
     */
    public final void setSigma(double s) {
        sigma = s;
        sigmaSquared = s*s;
    }
    public Dimension getSigmaDimension() {return Length.DIMENSION;}

    /**
     * Accessor method for Lennard-Jones energy parameter
     */
    public double getEpsilon() {return epsilon;}
    /**
     * Mutator method for Lennard-Jones energy parameter
     */
    public final void setEpsilon(double eps) {
        epsilon = eps;
        epsilon4 = eps*4.0;
        epsilon48 = eps*48.0;
        epsilon624 = eps*624.0;
    }
    public Dimension getEpsilonDimension() {return Energy.DIMENSION;}

    private double sigma, sigmaSquared;
    private double epsilon;
    private double epsilon4, epsilon48, epsilon624;
    private static final double _168div624 = 168./624.;

}
