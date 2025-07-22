package etomica.potential.COMPASS;


import etomica.potential.IPotentialAngleAngleHybrid;

public class AngleAngleCOMPASS implements IPotentialAngleAngleHybrid {
    private double k1;// Spring constant gives a measure of the strength of harmonic interaction

    private double theta20, theta10;

    public AngleAngleCOMPASS(double theta10, double theta20, double k1t){
        this.theta20 = theta20;
        this.theta10 = theta10;
        this.k1 = k1t;
    }

    public void u012add(double costheta1, double costheta2, double[] u012){
        double theta1 = Math.acos(costheta1);
        double theta2 = Math.acos(costheta2);
        u012[0] = k1 * (theta2-theta20) * (theta1-theta10);
    }

    public double u(double costheta1, double costheta2) {
        double theta1 = Math.acos(costheta1);
        double theta2 = Math.acos(costheta2);
        return  k1 * (theta2-theta20) * (theta1-theta10);
    }
}
