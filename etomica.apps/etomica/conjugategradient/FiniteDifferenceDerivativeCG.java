package etomica.conjugategradient;

import etomica.action.Activity;
import etomica.atom.AtomAgentManager;
import etomica.atom.IAtom;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.paracetamol.AnalyticalDerivativeEnergyParacetamol;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.util.numerical.FiniteDifferenceDerivative;

public class FiniteDifferenceDerivativeCG {
	
	/*
	 * This class is developed to calculate for second derivative of energy
	 *  by applying finite difference derivative to derivative function,
	 *  which the first derivative is determined by calculating:
	 *  	a. the force of molecule (1st derivative for translation)
	 *      b. the torque of molecule (1st derivative for rotation)
	 * 
	 * @author Tai Tan
	 */
	
	protected Box box;
	protected MeterPotentialEnergy meterEnergy;
	protected PotentialMaster potentialMaster;
	protected IteratorDirective allAtoms;
	protected PotentialCalculationForceSum forceSum;
	protected AtomAgentManager agentManager;
	protected Activity activity;
	
	protected AnalyticalDerivativeEnergyParacetamol derivativeFunction;
	protected double h;
	protected boolean hOptimizer;
	
	public FiniteDifferenceDerivativeCG(Box box, PotentialMaster potentialMaster, AnalyticalDerivativeEnergyParacetamol derivativeFunction){
		this.box = box;
		this.potentialMaster = potentialMaster;
		this.derivativeFunction = derivativeFunction;
		meterEnergy = new MeterPotentialEnergy(potentialMaster);
		allAtoms = new IteratorDirective();
		forceSum = new PotentialCalculationForceSum();
		h = 0.00001;
		hOptimizer = false;
		
		MyAgentSource source = new MyAgentSource(box.getSpace());
		agentManager = new AtomAgentManager(source, box);
		forceSum.setAgentManager(agentManager);
	}
	
	
	public double f(double[] u){
		return derivativeFunction.f(u);
	}
	
	public double df(int[] d, double[] u){
		return derivativeFunction.df(d, u);
	}
	
	public double[] d2f(int[] d, double[] u){
		
		int coordinateDim = u.length;
		double[] d2fdu2 = new double[coordinateDim];
		
		derivativeFunction = new AnalyticalDerivativeEnergyParacetamol(box, potentialMaster);
		double[] uDerivative = new double[coordinateDim];
		int[] dAssign = new int[coordinateDim];
		
		/*
		 *  To compute the first-order derivative from DerivativeFunctionParacetamol
		 *   for which uDerivative will be next input as u to compute for the second derivative
		 */
		for (int i=0; i< coordinateDim; i++){
			
			/*
			 *  To assign d to differentiate over all the dimensions
			 */
			for (int j=0; j< coordinateDim; j++){
				if (i==j){
					dAssign[j] = 1;
				} else{
					dAssign[j] = 0;
				}
			}
			uDerivative[i] = derivativeFunction.df(dAssign, u);
		}
		
		FiniteDifferenceDerivative finiteDifferenceDerivative = new FiniteDifferenceDerivative(derivativeFunction);
		finiteDifferenceDerivative.setH(0.000001);
		finiteDifferenceDerivative.setHOptimizer(true);
			
		/*
		 * Using Finite Difference Derivative to compute for the second derivative
		 *  of the energy.
		 */
			
		for (int i=0; i< coordinateDim; i++){
			
			/*
			 *  To assign d to differentiate over all the dimensions
			 */
			for (int j=0; j< coordinateDim; j++){
				if (i==j){
					dAssign[j] = 1;
				} else{
					dAssign[j] = 0;
				}
			}
			d2fdu2[i] = finiteDifferenceDerivative.df(dAssign, uDerivative);
		}
		
		return d2fdu2;
	}

	
	public static class MyAgentSource implements AgentSource{
		public MyAgentSource(Space space){
			this.space = space;
		}
		public void releaseAgent(Object agent, IAtom atom){}
		
		public Class getAgentClass(){
			return IntegratorVelocityVerlet.MyAgent.class;
		}
		
		public Object makeAgent(IAtom atom){
			return new IntegratorVelocityVerlet.MyAgent(space);
		}
		protected Space space;
	}


	public double getH() {
		return h;
	}

	public void setH(double h) {
		this.h = h;
	}

	public boolean isHOptimizer() {
		return hOptimizer;
	}

	public void setHOptimizer(boolean optimizer) {
		hOptimizer = optimizer;
	}
}
