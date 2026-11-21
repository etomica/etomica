package etomica.mappedDensity;

import etomica.atom.IAtom;
import etomica.potential.IPotential1;
import etomica.space.Vector;

/**
 *
 * parabolic external field
 */
public class P1Parabolic implements IPotential1 {
    double arg = 1;  //Vo

    @Override
    public double u(IAtom atom) {
        double z = atom.getPosition().getX(2);
        return arg*(z*z);
    }

    @Override
    public double udu(IAtom atom, Vector f) {
        double z = atom.getPosition().getX(2);
        f.setX(2, f.getX(2) + arg*2*(z));
        return arg*(z*z);
    }
}
