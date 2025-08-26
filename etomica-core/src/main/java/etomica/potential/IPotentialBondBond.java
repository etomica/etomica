package etomica.potential;

import etomica.atom.IAtom;
import etomica.exception.MethodNotImplementedException;
import etomica.space.Vector;

public interface IPotentialBondBond {
    default double u(double r10, double r11) {
        throw new MethodNotImplementedException();
    }


    default double du(double r10, double r11 ) {
        throw new MethodNotImplementedException();
    }


    default void u012add(double r10, double r11, double[] u012) {
        u012[0] += u(r10, r11);
        u012[1] += du(r10, r11);
        u012[2] += d2u(r10, r11);
    }


    default double d2u(double r10, double r11){
        throw new MethodNotImplementedException();
    }


    default double u(Vector dr12, Vector dr23, IAtom atom1, IAtom atom2, IAtom atom3)
    {
        return u(Math.sqrt(dr12.squared()),  Math.sqrt(dr23.squared()));
    }
}
