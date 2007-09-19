package etomica.virial;

import etomica.atom.AtomSet;
import etomica.box.Box;
import etomica.potential.IPotential;

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
	public double f(AtomSet pair, double beta);

	/**
	 * @return
	 */
	public IPotential getPotential();
	
	public void setBox(Box box);
}
