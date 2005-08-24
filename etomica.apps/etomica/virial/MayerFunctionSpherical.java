package etomica.virial;

import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.space.CoordinatePair;
import etomica.space.Space;

/**
 * @author kofke
 *
 * Abstract class for a Mayer f-function, which takes a pair of atoms and
 * returns exp(-u(pair)/kT) - 1
 */
public abstract class MayerFunctionSpherical implements MayerFunction, java.io.Serializable {

    public MayerFunctionSpherical(Space space) {
        coordPair = new CoordinatePair(space);
    }

    /**
     * returns Mayer function between atoms in the pair at temperature
     * 1/beta
     */
    public abstract double f(CoordinatePair cPair, double beta);

    public double f(AtomSet pair, double beta) {
        coordPair.reset((AtomPair)pair);
        return f(coordPair,beta);
    }

    private final CoordinatePair coordPair;
}
