package etomica.virial.dos;

import etomica.virial.MayerFunction;

/**
 * @author kofke
 *
 * Pair weight based on the density-of-states for the LJ model.  Probability for
 * a configuration is the reciprocal of the value of the density of states
 * evaluated at the configuration's energy.
 */
public class MayerDOS2 extends MayerFunction {
	
	/**
	 * Constructor for A2.
	 */
	public MayerDOS2() {
		super();
	}

	/**
	 * @see etomica.virial.MayerFunction#f(etomica.AtomPair, double)
	 */
	public double f(AtomPair pair, double beta) {
		double r2 = pair.r2();
		double r6 = r2*r2*r2;
		double u = 4/r6*(1/r6 - 1);
		double e = Math.exp(-u);
		double pi = pr(r6M(u))/e + ((e < 1) ? 0 : pr(r6P(u))/e);
		return 1.0/pi;
	}
	
	private double r6M(double u0) {
		return 2*(-1 + Math.sqrt(1 + u0))/u0;
	}
	
	private double r6P(double u0) {
		return 2*(-1 - Math.sqrt(1 + u0))/u0;
	}
	
	private double pr(double r6) {
		double pr = 4*Math.PI*Math.sqrt(r6)/(24./r6*(2.0/r6 - 1.0));
		return (pr<0) ? -pr : +pr;
	}

}