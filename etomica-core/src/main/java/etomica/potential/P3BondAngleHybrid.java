package etomica.potential;

public interface P3BondAngleHybrid {
    double u (double costheta, double r10);

    default double u (double costheta, double r12, double r12Sq, double r23Sq){
        return u(costheta, r12);
    }

    void udu(double costheta, double r12, double[] u, double[] du);
}
