/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential.amoeba;

import etomica.potential.IPotentialBondAngle;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;

/**
 * Simple 3-body soft bond-angle potential
 *
 * @author andrew
 */
public class P3AmeobaBondAngle implements IPotentialBondAngle {

    protected double angle;
    protected double epsilon;
    protected double c12, c23;
    protected double bl12, bl23;

    public P3AmeobaBondAngle(double angle, double epsilon) {
        this(angle, epsilon, 0, 0, 0, 0);
    }

    public P3AmeobaBondAngle(double angle, double epsilon, double c12, double c23, double bl12, double bl23) {
        setAngle(angle);
        setEpsilon(epsilon);
        setStretchBend(c12, c23, bl12, bl23);
    }

    public double u(double costheta) {
        double theta;
        if (costheta > 1) {
            theta = 0;
        } else if (costheta < -1) {
            theta = Math.PI;
        } else {
            theta = Math.acos(costheta);
        }
        double dtheta = theta - angle;
        return epsilon*dtheta*dtheta*(1 + dtheta*(-0.014 + dtheta*
                (5.6e-5 + dtheta*(-7e-7 + 2.2e-8*dtheta))));
    }

    public double u(double costheta, double r12Sq, double r23Sq) {
        double theta;
        if (costheta > 1) {
            theta = 0;
        } else if (costheta < -1) {
            theta = Math.PI;
        } else {
            theta = Math.acos(costheta);
        }
        double dtheta = theta - angle;
        double u = epsilon*dtheta*dtheta*(1 + dtheta*(-0.014 + dtheta*
                (5.6e-5 + dtheta*(-7e-7 + 2.2e-8*dtheta))));
        if (c12 == 0 && c23 == 0) return u;
        u += (c12*(Math.sqrt(r12Sq)-bl12) + c23*(Math.sqrt(r23Sq)-bl23)) * dtheta;
        return u;
    }

    /**
     * Sets the nominal bond angle (in radians)
     */
    public void setAngle(double newAngle) {
        angle = newAngle;
    }
    
    /**
     * Returns the nominal bond angle (in radians)
     */
    public double getAngle() {
        return angle;
    }

    public void setStretchBend(double c12, double c23, double bl12, double bl23) {
        this.c12 = c12;
        this.c23 = c23;
        this.bl12 = bl12;
        this.bl23 = bl23;
    }

    /**
     * Sets the characteristic energy of the potential
     */
    public void setEpsilon(double newEpsilon) {
        epsilon = newEpsilon;
    }
    
    /**
     * Returns the characteristic energy of the potential
     */
    public double getEpsilon() {
        return epsilon;
    }
    
    public Dimension getEpsilonDimension() {
        return Energy.DIMENSION;
    }
    
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public void udu(double costheta, double[] u, double[] du) {
        double theta;
        if (costheta > 1) {
            theta = 0;
        } else if (costheta < -1) {
            theta = Math.PI;
        } else {
            theta = Math.acos(costheta);
        }
        double r2d = 180/Math.PI;
        double dtheta = (theta - angle)*r2d;
        double sintheta = Math.sqrt(1 - costheta * costheta);
        du[0] = -epsilon * dtheta * (2 + dtheta*(3*-0.014 + dtheta*
                (4*5.6e-5 + dtheta*(6*-7e-7 + 2.2e-8*dtheta)))) / sintheta;
        u[0] = epsilon*dtheta*dtheta/(r2d*r2d)*(1 + dtheta*(-0.014 + dtheta*
                     (5.6e-5 + dtheta*(-7e-7 + 2.2e-8*dtheta))));

//        Unit kcalpmole = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT), Mole.UNIT);
//        System.out.println("angle "+theta+" "+dtheta+" "+kcalpmole.fromSim(epsilon)+" "+kcalpmole.fromSim(u[0]));
    }
}
