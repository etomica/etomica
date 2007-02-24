package etomica.virial;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.space.Space;
import etomica.space.Vector;

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

    public double f(AtomSet pair, double beta) {
        dr.Ev1Mv2(((AtomLeaf)((AtomPair)pair).atom1).getCoord().getPosition(),((AtomLeaf)((AtomPair)pair).atom0).getCoord().getPosition());
        return f(dr.squared(), beta);
    }

    protected final Vector dr;
}
