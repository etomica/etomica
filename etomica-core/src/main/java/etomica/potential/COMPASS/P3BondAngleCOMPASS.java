package etomica.potential.COMPASS;

import etomica.potential.IPotential2;
import etomica.potential.IPotentialBondAngle;
import etomica.potential.P3BondAngle;
import etomica.potential.P3BondAngleHybrid;

public class P3BondAngleCOMPASS implements IPotentialBondAngle {
    private double k2 = 100.0, k3=100, k4 =100;// Spring constant gives a measure of the strength of harmonic interaction
    public double theta0;

    public P3BondAngleCOMPASS(double k2, double k3, double k4, double theta0) {
        this.k2 = k2;
        this.k3 = k3;
        this.k4 = k4;
        this.theta0 = theta0;
    }

    public double u(double costheta){
        double theta= Math.acos(costheta);
        double dt = theta - theta0;
        return k2 * dt * dt + k3 * dt * dt * dt + k4 * dt * dt * dt * dt;
    }

    public void udu(double costheta, double[] u, double[] du){
        double theta = Math.acos(costheta);
        double dt = theta - theta0;
        u[0] = k2 * dt * dt + k3 * dt * dt * dt + k4 * dt * dt * dt * dt;
        du[0] = -2*k2*dt -3*k3*dt*dt -4*k4*dt*dt*dt;
    }
}
