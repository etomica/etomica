/*
 * Created on Apr 22, 2003
 *
 */
package etomica.dpd;

import etomica.AtomPair;
import etomica.Bond;
import etomica.EtomicaElement;
import etomica.Potential2SoftSpherical;
import etomica.Simulation;
import etomica.SimulationElement;
import etomica.Space;
import etomica.Default;
import etomica.Potential2;

/**
 * @author cribbin
 *
 * The potential used for dissipative particle dynamic simulations.
 * Based on equation 1 of Willemsen et al. (J. Comp. Phys. 162 pp.385-394 (2000))
 */
public class P2DPD extends Potential2 implements Potential2.Soft, EtomicaElement {
	double maxRepel; 	//maximum repulsion parameter (aij)
	double sigma;		//noise amplitude parameter
	double rC, rC2;			//cutoff radius, square
	double gamma;
	double temperature = Default.TEMPERATURE;
	double timeStep = Default.TIME_STEP;
	double constR;      //temperature-dependent constant appearing in random force
	Space.Vector rHat, fDissipative, fConservative, fRandom, grad;
	
	public String getVersion() {return "P2DPD:03.04.25/"+Potential2SoftSpherical.VERSION;}
	
	public P2DPD(){
		this(3.0, 1.0, 1.0);
	}
	
	public P2DPD(SimulationElement parent) {
		this(parent, 3.0, 1.0, 1.0);
	}

	public P2DPD(double sigma, double maxRepel, double cutoffRadius){
		this(Simulation.instance.hamiltonian.potential, sigma, maxRepel, cutoffRadius);
	}

	public P2DPD(SimulationElement parent, double sigma, double maxRepel, double cutoffRadius){
		super(parent);
		setSigma(sigma);
		setMaxRepel(maxRepel);
		setRC(cutoffRadius);
		rHat = space.makeVector();
		fDissipative = space.makeVector();
		fConservative = space.makeVector();
		fRandom = space.makeVector();
		grad = space.makeVector();
	}

	/**
	 * Not implemented
	 * @see etomica.Potential2#energy(AtomPair)
	 */
	public double energy(AtomPair pair) {
		throw new etomica.exception.MethodNotImplementedException();
	}
	/**
	 * Not implemented
	 * @see etomica.Potential2.Soft#virial(AtomPair)
	 */
	public double virial(AtomPair pair) {
		throw new etomica.exception.MethodNotImplementedException();
	}
	/**
	 * Not implemented
	 * @see etomica.Potential2.Soft#hyperVirial(AtomPair)
	 */
	public double hyperVirial(AtomPair pair) {
		throw new etomica.exception.MethodNotImplementedException();
	}
	/**
	 * Not implemented
	 * @see etomica.Potential2.Soft#integral(double)
	 */
	public double integral(double rCut) {
		throw new etomica.exception.MethodNotImplementedException();
	}

	/**
	 * Returns the gradient of the potential. 
	 * Overwrites base class method.
	 * @param an AtomPair; the two atoms the gradient is calculated for.
	 * @return The gradient of the vector between a pair of atoms.
	 */
	public Space.Vector gradient(AtomPair pair) {
		double r2 = pair.r2();
		if(r2 < rC2) {
			pair.cPair.resetV();
			double r = Math.sqrt(r2);
			double rand = Math.sqrt(3.0)*(2.*Simulation.random.nextDouble()-1.);
//			double rand = (2.*Simulation.random.nextDouble()-1.);
//			double rand = Simulation.random.nextGaussian();

			double wR = 1.0 - r/rC;
			double g = 0.0;
//			g += maxRepel(pair);  //conservative     maxRepel*cR/r dr
//			g += -gamma*wR*pair.vDotr()/r;  //dissipative  -gamma*cR*cR*vdotr/r dr
//			g += constR*rand;// random  constR*cR/r dr,  constR = sigma/sqrt(dt)
			grad.Ea1Tv1(-g*wR/r, pair.dr());  //minus because force is negative of gradient
//			grad.TE(-1.0);
//			fConservative.Ea1Tv1(maxRepel*cR/r, pair.dr());
//			fDissipative.Ea1Tv1(-gamma*cR*cR*pair.vDotr()/r,pair.dr());
//			if(Math.abs(rand)>0.5) System.out.println(rand);			
//			fRandom.Ea1Tv1(constR*cR/r*rand,pair.dr());
//			grad.E(0.0);
//			grad.E(fConservative);
//			grad.PE(fDissipative);
//			grad.PE(fRandom);
		} else {
			grad.E(0.0);
		} 	
		double r0 = 0.5;
		if(Bond.areBonded(pair.atom1, pair.atom2))	{
			grad.PEa1Tv1(5.*(1.0 - r0/Math.sqrt(pair.r2())),pair.dr());
		}
		return grad;
	}
	
	protected double maxRepel(AtomPair pair) {
		return maxRepel;
	}

	/**
	 * @return the maximum repulsion parameter
	 */
	public double getMaxRepel() {
		return maxRepel;
	}

	/**
	 * @return the cutoff radius
	 */
	public double getRC() {
		return rC;
	}

	/**
	 * @return the noise amplitude
	 */
	public double getSigma() {
		return sigma;
	}

	/**
	 * sets the maximum repulsion parameter
	 * @param d
	 */
	public void setMaxRepel(double d) {
		maxRepel = d;
	}

	/**
	 * sets the cutoff radius
	 * @param d
	 */
	public void setRC(double d) {
		rC = d;
		rC2 = d*d;
	}

	/**
	 * sets the noise amplitude
	 * @param d
	 */
	public void setSigma(double d) {
		sigma = d;
		updateConstants();
	}
	
	public void setTemperature(double t) {
		temperature = t;
		updateConstants();
	}
	
	private void updateConstants() {
		gamma = 0.5*sigma*sigma/temperature;
		constR = sigma/Math.sqrt(timeStep);
	}



	/**
	 * Returns the timeStep.
	 * @return double
	 */
	public double getTimeStep() {
		return timeStep;
	}

	/**
	 * Sets the timeStep.
	 * @param timeStep The timeStep to set
	 */
	public void setTimeStep(double timeStep) {
		this.timeStep = timeStep;
		updateConstants();
	}

}//end class-P2DPD