package etomica.virial;

import etomica.api.IBox;
import etomica.api.IMoleculeList;

/**
 * Interface for a Mayer f-function, which takes 3+ atoms and
 * returns exp(-(u(atoms)-u(pairs))/kT) - 1
 *
 * @author Andrew Schultz
 */
public interface MayerFunctionNonAdditive {

    /**
     * returns exp(-beta*(U - Upair))
     * 
     * r2 is listed in the order
     *   (0,1),(0,2)...(0,n-1),(1,2),(1,3)...(1,n-1)...(n-2,n-1)
     */
	public double f(IMoleculeList molecules, double[] r2, double beta);

	public void setBox(IBox box);
}
