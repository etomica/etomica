package etomica.virial;

import etomica.AtomPair;

/**
 * @author kofke
 *
 * Abstract class for a Mayer f-function, which takes a pair of atoms and
 * returns exp(-u(pair)/kT) - 1
 */
public abstract class MayerFunction {

	/**
	 * Constructor for MayerFunction.
	 */
	public MayerFunction() {
		super();
	}

	public abstract double f(AtomPair pair, double beta);
	
	public String toString() {return "f   ";}
}
