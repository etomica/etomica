package etomica.potential;

import etomica.atom.IAtom;
import etomica.exception.MethodNotImplementedException;
import etomica.space.Vector;

public interface IPotentialBondTorsionHybrid {
    default double u(double r1, double costheta) {
        throw new MethodNotImplementedException();
    }

    default double du(double r1, double costheta) {
        throw new MethodNotImplementedException();
    }

    default double d2u(double r1, double costheta) {
        throw new MethodNotImplementedException();
    }

    default void u012add(double r, double costheta, double[] u012) {
        u012[0] += u(r, costheta);
        u012[1] += du(r, costheta);
        u012[2] += d2u(r, costheta);
    }

    default double u(Vector dr12, Vector dr23, double costheta, IAtom atom1, IAtom atom2, IAtom atom3){
        dr12.Ev1Mv2(dr12, dr23);
        return u(Math.sqrt(dr12.squared()), costheta);
    }

}
