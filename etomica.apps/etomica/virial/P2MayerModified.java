package etomica.virial;

import etomica.AtomPair;
import etomica.Potential2;
import etomica.PotentialTruncation;
import etomica.SimulationElement;

/**
 * @author David Kofke
 *
 * Mayer "potential", modified to add a core hard sphere.
 */
public class P2MayerModified extends Potential2 {

	/**
	 * Constructor for P2MayerModified.
	 * @param parent
	 */
	public P2MayerModified(SimulationElement parent) {
		super(parent);
	}

	/**
	 * Constructor for P2MayerModified.
	 * @param parent
	 * @param potentialTruncation
	 */
	public P2MayerModified(
		SimulationElement parent,
		PotentialTruncation potentialTruncation) {
		super(parent, potentialTruncation);
	}

	/**
	 * Returns -T*ln(abs(f)), where f is the Mayer function exp(-u/T)-1, or if r
	 * is less than core hard-sphere diameter, return max of -T*ln(abs(f)) and
	 * 1.0.
	 * @see etomica.Potential2#energy(etomica.AtomPair)
	 */
	public double energy(AtomPair pair) {
		double r2 = pair.r2();
		bu = beta*potential.energy(pair);
		if(r2 < sigma2 && bu > UF1) return 0.0;
		
		double f = Math.exp(-bu) - 1.0;		
		return -temperature*Math.log((f>0)?f:-f);  //argument to log is abs(f)
	}
	
	public double mostRecentBetaU() {
		return bu;
	}

	/**
	 * Returns sigma, the diameter of the core hard sphere.
	 * @return double
	 */
	public double getSigma() {
		return sigma;
	}

	/**
	 * Returns the temperature.
	 * @return double
	 */
	public double getTemperature() {
		return temperature;
	}

	/**
	 * Sets sigma, the core hard-sphere diameter.
	 * @param sigma The sigma to set
	 */
	public void setSigma(double sigma) {
		this.sigma = sigma;
		sigma2 = sigma*sigma;
	}

	/**
	 * Sets the temperature.
	 * @param temperature The temperature to set
	 */
	public void setTemperature(double temperature) {
		this.temperature = temperature;
		beta = 1.0/temperature;
	}

	/**
	 * Returns the potential.
	 * @return Potential2
	 */
	public Potential2 getPotential() {
		return potential;
	}

	/**
	 * Sets the potential.
	 * @param potential The potential to set
	 */
	public void setPotential(Potential2 potential) {
		this.potential = potential;
	}

	private double temperature, beta, sigma, sigma2, bu;
	private Potential2 potential;
	private static final double UF1 = -Math.log(2.0); //value of beta*u for which f = +1

}
