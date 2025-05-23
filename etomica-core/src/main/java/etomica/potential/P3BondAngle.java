/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Boundary;
import etomica.units.Degree;
import etomica.units.dimensions.Angle;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;

/**
 * Simple 3-body soft bond-angle potential
 *
 * @author andrew
 */
public class P3BondAngle implements IPotentialBondAngle {

    public P3BondAngle(double angle, double epsilon) {
        setAngle(angle);
        setEpsilon(epsilon);
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
        return 0.5*epsilon*dtheta*dtheta;
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
    
    public Dimension getAngleDimension() {
        return Angle.DIMENSION;
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
        double dtheta = theta - angle;
//        if (Math.abs(Degree.UNIT.fromSim(angle) - 112) < 0.01) {
//            System.out.println(Degree.UNIT.fromSim(theta));
//        }
        du[0] = -epsilon * dtheta / Math.sqrt(1 - costheta * costheta);
        u[0] = 0.5 * epsilon * dtheta * dtheta;
    }

    protected Boundary boundary;
    protected double angle;
    protected double epsilon;
}
