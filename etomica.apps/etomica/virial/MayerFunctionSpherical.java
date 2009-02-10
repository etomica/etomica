package etomica.virial;

import etomica.api.IMoleculeList;
import etomica.api.IVectorMutable;
import etomica.space.ISpace;

/**
 * @author kofke
 *
 * Abstract class for a Mayer f-function, which takes a pair of atoms and
 * returns exp(-u(pair)/kT) - 1
 */
public abstract class MayerFunctionSpherical implements MayerFunction, java.io.Serializable {

    public MayerFunctionSpherical(ISpace space) {
        dr = space.makeVector();
    }

    /**
     * returns Mayer function between atoms in the pair at temperature
     * 1/beta
     */
    public abstract double f(double r2, double beta);

    public double f(IMoleculeList pair, double r2, double beta) {
        return f(r2, beta);
    }
    
    protected final IVectorMutable dr;
}
