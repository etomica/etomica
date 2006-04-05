package etomica.virial;

import etomica.simulation.Simulation;
import etomica.space.CoordinatePair;
import etomica.space.Space;

/**
 * @author kofke
 *
 * Hard-sphere Mayer function.  -1 if r < sigma; 0 otherwise
 */
public class MayerHardSphere extends MayerFunctionSpherical {

	private double sigma, sigma2;
	/**
	 * Constructor for MayerHardSphere.
	 */
	public MayerHardSphere(Simulation sim) {
		this(sim.space, sim.getDefaults().atomSize);
	}
	public MayerHardSphere(Space space, double sigma) {
        super(space);
		setSigma(sigma);
	}

	/**
	 * @see etomica.virial.MayerFunctionSpherical#f(etomica.AtomPair)
	 */
	public double f(double r2, double beta) {
		return (r2<sigma2) ? -1.0 : 0.0;
	}

	/**
	 * Returns the HS diameter.
	 * @return double
	 */
	public double getSigma() {
		return sigma;
	}

	/**
	 * Sets the HS diameter.
	 * @param sigma The sigma to set
	 */
	public void setSigma(double sigma) {
		this.sigma = sigma;
		sigma2 = sigma*sigma;
	}

}
