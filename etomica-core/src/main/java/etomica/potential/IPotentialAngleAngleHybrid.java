package etomica.potential;

import etomica.atom.IAtom;
import etomica.exception.MethodNotImplementedException;
import etomica.space.Vector;

public interface IPotentialAngleAngleHybrid {
    default double u(double r1, double costheta2) {
        throw new MethodNotImplementedException();
    }


    default double du(double costheta1, double costheta2) {
        throw new MethodNotImplementedException();
    }


    default double d2u(double costheta1, double costheta2) {
        throw new MethodNotImplementedException();
    }


    default void u012add(double costheta1, double costheta2, double[] u012) {
        u012[0] += u(costheta1, costheta2);
        u012[1] += du(costheta1, costheta2);
        u012[2] += d2u(costheta1, costheta2);
    }


    default double u(Vector dr12, Vector dr23, double costheta, IAtom atom1, IAtom atom2, IAtom atom3){
        dr12.Ev1Mv2(dr12, dr23);
        return u(Math.sqrt(dr12.squared()), costheta);
    }


}


