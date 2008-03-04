package etomica.virial;

import etomica.api.IVector;
import etomica.api.IAtomSet;
import etomica.api.IAtomPositioned;

import etomica.space.Space;

/**
 * @author kofke
 *
 * Abstract class for a Mayer f-function, which takes a pair of atoms and
 * returns exp(-u(pair)/kT) - 1
 */
public abstract class MayerFunctionSpherical implements MayerFunction, java.io.Serializable {

    public MayerFunctionSpherical(Space space) {
        dr = space.makeVector();
    }

    /**
     * returns Mayer function between atoms in the pair at temperature
     * 1/beta
     */
    public abstract double f(double r2, double beta);

    public double f(IAtomSet pair, double beta) {
        dr.Ev1Mv2(((IAtomPositioned)pair.getAtom(1)).getPosition(),((IAtomPositioned)pair.getAtom(0)).getPosition());
        return f(dr.squared(), beta);
    }
    
    protected final IVector dr;
}
