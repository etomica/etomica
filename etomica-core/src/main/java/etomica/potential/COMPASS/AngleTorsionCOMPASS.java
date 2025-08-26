package etomica.potential.COMPASS;

import etomica.potential.IPotentialAngleTorsionHybrid;

public class AngleTorsionCOMPASS implements IPotentialAngleTorsionHybrid {
    private double k1, k2, k3;// Spring constant gives a measure of the strength of harmonic interaction

    private double phi0, theta0;

    public AngleTorsionCOMPASS(double theta0, double phi0, double k1t, double k2t, double k3t){
        this.theta0 = theta0;
        this.phi0 = phi0;
        this.k1 = k1t;
        this.k2 = k2t;
        this.k3 = k3t;
    }

    public void u012add(double costheta, double cosphi, double[] u012){
        double phi = Math.acos(cosphi);
        double theta = Math.acos(costheta);
        double torsion = (k1 * cosphi) + (k2 * Math.cos(2*phi)) + (k3 * Math.cos(3*phi));
        u012[0] = (theta-theta0) * torsion;
    }

    public double u(double costheta, double cosphi) {
        double phi = Math.acos(cosphi);
        double theta = Math.acos(costheta);
        double torsion = (k1 * cosphi) + (k2 * Math.cos(2*phi)) + (k3 * Math.cos(3*phi));
        return (theta-theta0) * torsion;
    }


}
