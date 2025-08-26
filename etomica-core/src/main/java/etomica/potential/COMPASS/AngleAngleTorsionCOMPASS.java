package etomica.potential.COMPASS;

import etomica.potential.IPotentialAngleAngleTorsionHybrid;

public class AngleAngleTorsionCOMPASS implements IPotentialAngleAngleTorsionHybrid {
    private double k1;// Spring constant gives a measure of the strength of harmonic interaction

    private double phi0, theta0, theta10;

    public AngleAngleTorsionCOMPASS(double theta0, double theta10, double phi0, double k1t){
        this.theta0 = theta0;
        this.theta10 = theta10;
        this.phi0 = phi0;
        this.k1 = k1t;
    }

    public void u012add(double costheta, double costheta1,  double cosphi, double[] u012){
        double theta = Math.acos(costheta);
        double theta1 = Math.acos(costheta1);
        u012[0] = k1 * (theta-theta0) * (theta1-theta10) * cosphi;
    }

    public double u(double costheta, double costheta1,  double cosphi) {
        double theta = Math.acos(costheta);
        double theta1 = Math.acos(costheta1);
        return  k1 * (theta-theta0) * (theta1-theta10) * cosphi;
    }


}
