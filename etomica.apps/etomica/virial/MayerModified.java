package etomica.virial;

import etomica.AtomPair;
import etomica.Default;
import etomica.potential.Potential2;

/**
 * @author kofke
 *
 * A Mayer function that returns abs(f) or, when separation is less than a
 * specified amount (sigma), max (abs (f), 1.0).  Suitable for MC sampling.
 */
public class MayerModified extends MayerFunction {

	/**
	 * Constructor for MayerModified.
	 */
	public MayerModified(Potential2 potential) {
		this(potential, Default.ATOM_SIZE);
	}
	public MayerModified(Potential2 potential, double sigma) {
		super();
		this.potential = potential;
		setSigma(sigma);
	}

	/**
	 * @see etomica.virial.MayerFunction#f(etomica.AtomPair, double)
	 */
	public double f(AtomPair pair, double beta) {
		double r2 = pair.r2();
		double bu = beta*potential.energy(pair);
		if(r2 < sigma2 && bu > UF1) return 1.0;
//		if(r2 < sigma2 ) return 1.0;
//		double bu = beta*potential.energy(pair);
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

	private final Potential2 potential;
	private static final double UF1 = -Math.log(2.0); //value of beta*u for which f = +1
	private double sigma, sigma2;

}
