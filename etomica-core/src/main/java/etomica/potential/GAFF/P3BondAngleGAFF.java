package etomica.potential.GAFF;

import etomica.potential.IPotentialBondAngle;
import etomica.potential.P3BondAngle;
import etomica.units.*;

public class P3BondAngleGAFF implements IPotentialBondAngle {
    public P3BondAngleGAFF(double epsilon, double angle) {
        setAngle(angle);
        setEpsilon(epsilon);
    }

    public void setAngle(double newAngle) {
        angle = newAngle;
    }

    public void setEpsilon(double newEpsilon) {
        epsilon = newEpsilon;
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
        return epsilon*dtheta*dtheta;
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
        Unit kjoulepmole = new UnitRatio(new PrefixedUnit(Prefix.KILO, Joule.UNIT), Mole.UNIT);
        double wnew = kjoulepmole.fromSim(epsilon);
        double dx = theta - angle;
        double denergy = wnew*dx*dx;
        System.out.println(dx + " " + theta +" " + angle + " " + wnew + " Energy : " + denergy );
        u[0] =  epsilon * dtheta * dtheta;
        du[0] = 2 * epsilon * dtheta;
    }
    public double angle, epsilon;
}
