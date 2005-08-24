package etomica.virial;

import etomica.space.CoordinatePair;
import etomica.space.Space;
import etomica.util.Default;

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
	public MayerHardSphere(Space space) {
		this(space, Default.ATOM_SIZE);
	}
	public MayerHardSphere(Space space, double sigma) {
        super(space);
		setSigma(sigma);
	}

	/**
	 * @see etomica.virial.MayerFunctionSpherical#f(etomica.AtomPair)
	 */
	public double f(CoordinatePair cPair, double beta) {
		return (cPair.r2()<sigma2) ? -1.0 : 0.0;
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
