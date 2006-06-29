package etomica.yukawa;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.potential.Potential2HardSpherical;
import etomica.simulation.Simulation;
import etomica.space.Space;

/**
 * Hard-core plus two Yukawa fluid (HC2Yukawa): A Lennard Jones like potential.
 * 
 * ============================================================================================================================
 * 2001. Pini, Stell, and Wilding. "Liquid-gas phase behaviour of an argon-like fluid modelled by the hard-core two-Yukawa potential"
 * ============================================================================================================================
 * 	
 * 			| infinity																		r <= sigma
 * U(r) =	| 
 * 			|(A1/r) * exp[-z1 * (r - sigma)] - (A2 * epsilon/r) * exp[-z2 * (r - sigma)]	r > sigma
 * 
 * 	LJ behavior parameters:	A1 = 1.6438 * sigma
 * 							z1 = 14.7 / sigma
 * 							A2 = 2.03 * sigma
 * 							z2 = 2.69 / sigma
 * 
 * @author msellers
 */

public final class P2HC2Yukawa extends Potential2HardSpherical implements EtomicaElement {

	public P2HC2Yukawa(Simulation sim){
		this(sim.space, sim.getDefaults().atomSize, sim.getDefaults().potentialWell);
	}
	
	public P2HC2Yukawa(Space space, double sigma, double epsilon){
		super(space);
		setSigma(sigma);
		setEpsilon(epsilon);
		setParameters(sigma);
	}
	
	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo("Hard-core plus two Yukawa fluid potential");
	}
	
	/**
	 * Energy method.  u(double r) form.
	 */
	public double u(double r2){
		double r = Math.sqrt(r2);
		//hard core repulsion
		if (r <= sigma){
			return Double.POSITIVE_INFINITY;
		}
		//two-Yukawa tail attraction
		if (r > sigma){
			
			if (r != rLast){
				term1 = (A1 / r) * Math.exp(-z1 * (r - sigma));
				term2 = (A2 * epsilon / r) * Math.exp(-z2 * (r - sigma));
				rLast = r;
			}
			
		return term1 - term2;
		}
	}
	
	/**
	 * First Derivative method. du/dr.
	 */
	public double du(double r2){
		double r = Math.sqrt(r2);
		//hard core repulsion
		if (r <= sigma){
			return Double.POSITIVE_INFINITY;
		}
		//two-Yukawa tail attraction
		if (r > sigma){
			
			if(r != rLast){
				term1 = (A1 / r) * Math.exp(-z1 * (r - sigma));
				term2 = (A2 * epsilon / r) * Math.exp(-z2 * (r - sigma));
				rLast = r;
			}
		return (-z1 * term1) + (z2 * term2);
		}
	}
	
	
	
	
	
	

	/**
	 * Accessor methods for size and energy parameters.
	 */
	public double getSigma() {return sigma;}
	
	public double getEpsilon() {return epsilon;}
	
	/**
	 * Mutator methods for size and energy parameters.
	 */
	public final void setEpsilon(double eps) {epsilon = eps;}
	
	public final void setSigma(double s) {sigma = s;}
	
	/**
	 * Parameter calculation method.
	 */
	public final void setParameters(double s){
		A1 = 1.6438 * s;
		A2 = 2.03 * s;
		z1 = 14.7 / s;
		z2 = 2.69 / s;
	}
	
	private double sigma;
	private double epsilon;
	private double A1;
	private double A2;
	private double z1;
	private double z2;
	private double rLast = -1.0;
	private double term1;
	private double term2;
	
}	
