package etomica.virial;

import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IPotential;

/**
 * @author kofke
 *
 * General Mayer function, which wraps the Mayer potential around an instance of
 * a Potential2 object.
 */
public class MayerEGeneral implements MayerFunction, java.io.Serializable {

	/**
	 * Constructor Mayer function using given potential.
	 */
	public MayerEGeneral(IPotential potential) {
		this.potential = potential;
	}

	public double f(IAtomSet pair, double beta) {
		return Math.exp(-beta*potential.energy(pair));
	}

	private final IPotential potential;

	/* (non-Javadoc)
	 * @see etomica.virial.MayerFunction#getPotential()
	 */
	public IPotential getPotential() {
		// TODO Auto-generated method stub
		return potential;
	}
	
	public void setBox(IBox newBox) {
	    potential.setBox(newBox);
	}
}
