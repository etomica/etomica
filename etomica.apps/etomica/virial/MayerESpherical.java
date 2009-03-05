package etomica.virial;

import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.IPotential;
import etomica.potential.Potential2Spherical;
/**
 * @author kofke
 *
 * Simple e-bond, exp(-beta*u).
 */
public class MayerESpherical implements MayerFunction {

	/**
	 * Constructor for MayerESpherical.
	 */
	public MayerESpherical(Potential2Spherical potential) {
		this.potential = potential;
	}

	/**
	 * @see etomica.virial.MayerFunctionSpherical#f(etomica.AtomPair, double, double)
	 */
	public double f(IMoleculeList pair, double r2, double beta) {
		return Math.exp(-beta*potential.u(r2));
	}
	
	public void setBox(IBox newBox) {
	    potential.setBox(newBox);
	}

	private final Potential2Spherical potential;

	public IPotential getPotential() {
		return potential;
	}

}
