package etomica.virial;

import etomica.potential.Potential2Spherical;
import etomica.simulation.Simulation;
import etomica.space.Space;

/**
 * @author kofke
 *
 * A Mayer function that returns abs(f) or, when separation is less than a
 * specified amount (sigma), max (abs (f), 1.0).  Suitable for MC sampling.
 */
public class MayerModified extends MayerFunctionSpherical {

	/**
	 * Constructor for MayerModified.
	 */
	public MayerModified(Simulation sim, Potential2Spherical potential) {
		this(sim.getSpace(), potential, sim.getDefaults().atomSize);
	}
	public MayerModified(Space space, Potential2Spherical potential, double sigma) {
        super(space);
		this.potential = potential;
		setSigma(sigma);
	}

	/**
	 * @see etomica.virial.MayerFunctionSpherical#f(etomica.AtomPair, double)
	 */
	public double f(double r2, double beta) {
		double bu = beta*potential.u(r2);
		if(r2 < sigma2 && bu > UF1) return 1.0;
		if(bu > -1e-6) return -bu*(1-0.5*bu); //repulsive region already eliminated, so bu close to zero if this is true
		double f = Math.exp(-bu) - 1.0;		
		return (f>0) ? f : -f;  
	}

	/**
	 * Returns sigma, the diameter of the core hard sphere.
	 * @return double
	 */
	public double getSigma() {
		return sigma;
	}

	/**
	 * Sets sigma, the core hard-sphere diameter.
	 * @param sigma The sigma to set
	 */
	public void setSigma(double sigma) {
		this.sigma = sigma;
		sigma2 = sigma*sigma;
	}

    private static final long serialVersionUID = 1L;
	private final Potential2Spherical potential;
	private static final double UF1 = -Math.log(2.0); //value of beta*u for which f = +1
	private double sigma, sigma2;

}
