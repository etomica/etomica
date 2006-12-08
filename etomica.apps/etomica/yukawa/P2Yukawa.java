package etomica.yukawa;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.potential.Potential2SoftSpherical;
import etomica.simulation.Simulation;
import etomica.space.Space;

/**
 * Yukawa interatomic potential.
 * 
 * U(r) = Vo * exp(-K * r) / r
 * 
 * where Vo is the potential's energy parameter and K is the characteristic distance (1/K is a measure for the screening length).
 * 
 * @author msellers
 *
 */

public final class P2Yukawa extends Potential2SoftSpherical implements EtomicaElement {
	
	public P2Yukawa(Simulation sim) {
		this(sim.getSpace(), sim.getDefaults().atomSize, sim.getDefaults().potentialWell);
	}
	
	public P2Yukawa(Space space, double kappa, double vzero){
		super(space);
		setKappa(kappa);
		setVZero(vzero);
	}
	
	public static EtomicaInfo getEtomicaInfo(){
		EtomicaInfo info = new EtomicaInfo("Yukawa potential.");
		return info;
	}
	
	/**
	 * Energy method.  u(double r) form.
	 */
	public double u(double r2){
		double r = Math.sqrt(r2);
		return vzero * Math.exp(-kappa * r) / r;
	}
	
	/**
	 * r * du/dr method.
	 */
	public double du(double r2){
		double r = Math.sqrt(r2);
		return (-vzero * Math.exp(-kappa * r)) * (kappa + (1 / r));
	}
	
	/**
	 * r^2 * d^2u/dr^2 method.
	 */
	public double d2u(double r2){
		double r = Math.sqrt(r2);
		return (vzero * Math.exp(-kappa * r) * r) * (kappa * (kappa + (2 / r)) + (2 / r2));
	}
	
	
    /**
     * Integral from rC to infinity.
     */
	public double uInt(double rC){
		
		return 0;
	}

	/**
	 * Accessor methods for size and energy parameters.
	 */
	public double getKappa() {return kappa;}
	
	public double getVZero() {return vzero;}
	
	/**
	 * Mutator methods for size and energy parameters.
	 */
	public final void setKappa(double k) {kappa = k;}
	
	public final void setVZero(double v) {vzero = v;}

	
	private double kappa;
	private double vzero;
}

