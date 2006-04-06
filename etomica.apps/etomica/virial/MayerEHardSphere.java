package etomica.virial;

import etomica.space.Space;

/**
 * @author kofke
 *
 * The e-function for hard spheres, returning 0 for r<sigma, 1 otherwise.
 */
public class MayerEHardSphere extends MayerESpherical {

	private double sigma, sigma2;

	/**
	 * Constructor for MayerEHardSphere.
	 * @param potential
	 */
	public MayerEHardSphere(Space space) {
		this(space,1.0);
	}
	
	public MayerEHardSphere(Space space, double sigma) {
		super(space,null);
		setSigma(sigma);
	}

	public double f(double r2, double beta) {
		return (r2<sigma2) ? 0.0 : 1.0;
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
