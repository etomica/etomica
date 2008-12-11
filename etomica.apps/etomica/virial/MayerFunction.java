package etomica.virial;

import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.IPotential;

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
	public double f(IMoleculeList pair, double beta);

	/**
	 * @return
	 */
	public IPotential getPotential();
	
	public void setBox(IBox box);
}
