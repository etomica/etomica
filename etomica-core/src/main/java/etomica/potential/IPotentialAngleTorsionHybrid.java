package etomica.potential;

import etomica.atom.IAtom;
import etomica.exception.MethodNotImplementedException;
import etomica.space.Vector;

public interface IPotentialAngleTorsionHybrid {
    default double u(double costheta, double cosphi) {
        throw new MethodNotImplementedException();
    }

    default double du(double costheta, double cosphi) {
        throw new MethodNotImplementedException();
    }

    default double d2u(double costheta, double cosphi) {
        throw new MethodNotImplementedException();
    }

    default void u012add(double costheta, double cosphi, double[] u012) {
        u012[0] += u(costheta, cosphi);
        u012[1] += du(costheta, cosphi);
        u012[2] += d2u(costheta, cosphi);
    }

    default double u(Vector dr12, Vector dr23, double costheta, IAtom atom1, IAtom atom2, IAtom atom3){
        dr12.Ev1Mv2(dr12, dr23);
        return u(Math.sqrt(dr12.squared()), costheta);
    }
}
