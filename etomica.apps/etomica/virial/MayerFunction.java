package etomica.virial;

import etomica.space.CoordinatePair;

/**
 * @author kofke
 *
 * Abstract class for a Mayer f-function, which takes a pair of atoms and
 * returns exp(-u(pair)/kT) - 1
 */
public interface MayerFunction {

    /**
     * returns Mayer function between atoms in the pair at temperature
     * 1/beta
     */
	public abstract double f(CoordinatePair cPair, double beta);
	
}
