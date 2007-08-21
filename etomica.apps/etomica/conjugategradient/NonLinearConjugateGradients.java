package etomica.conjugategradient;

import etomica.action.Activity;
import etomica.atom.AtomAgentManager;
import etomica.atom.IAtom;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntegratorHardField.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;

public class NonLinearConjugateGradients {

	/*
	 *   Nonlinear Conjugate Gradients with Newton-Raphson and Fletcher-Reeves
	 *  "An Introduction to the Conjugate Gradient Method w/o the Agonizing Pain"
	 *                      Jonathan Richard Shewchuk
	 *                           August 4, 1994
	 *                           
	 *  @author Tai Tan
	 */
	
	protected Box box;
	protected MeterPotentialEnergy meterEnergy;
	protected PotentialMaster potentialMaster;
	protected IteratorDirective allAtoms;
	protected PotentialCalculationForceSum forceSum;
	protected AtomAgentManager agentManager;
	protected Activity activity;
	protected FiniteDifferenceDerivative fDoublePrime;
	
	
	public NonLinearConjugateGradients(Box box, PotentialMaster potentialMaster){
		this.box = box;
		this.potentialMaster = potentialMaster;
		meterEnergy = new MeterPotentialEnergy(potentialMaster);
		allAtoms = new IteratorDirective();
		forceSum = new PotentialCalculationForceSum();
		
		MyAgentSource source = new MyAgentSource(box.getSpace());
		agentManager = new AtomAgentManager(source, box);
		forceSum.setAgentManager(agentManager);

	}
	
	public void NonLinearCG(fFunction function, int imax, int jmax, double epsilon, double[] u){
		
		/*
		 *  imax is a maximum number of CG iterations
		 *  jmax is a maximum number of Newton-Raphson iterations
		 *  
		 */
		
		int i=0;
		int k=0;
		int n;
		
		double deltaNew = 0;
		
		int coordinateDim = u.length;
		double[] fPrimeVal = function.fPrime(u);
		
		double[] r = new double[coordinateDim]; 
		double[] d = new double[coordinateDim];
		
		for (n=0; n<coordinateDim; n++){
			r[n] = - fPrimeVal[n];
			d[n] = r[n];
			
			System.out.println("r["+n + "] is: " + r[n]);
			deltaNew += r[n]*r[n];
		}
		
		double delta0 = deltaNew;
		double epsilon2_delta0 = epsilon*epsilon*delta0;
		
		System.out.println("NonlinearCG before WHILE loop...");
		System.out.println("imax is: "+ imax+" i is "+ i);
		System.out.println("DeltaNew: "+ deltaNew);
		System.out.println("epsilon: "+ epsilon);
		System.out.println("epsilon2_delta0: "+ epsilon2_delta0);
		
		while(i<imax && deltaNew > epsilon2_delta0){
			int j=0;
			double deltad = 0;
			double alpha_num = 0;
			double alpha_denom = 0;
			
			fPrimeVal = function.fPrime(u); 
			double[] fDoublePrimeVal = fDoublePrime.finiteDerivative(function, u, 0.001, coordinateDim, false);
			
			for(n=0; n<coordinateDim; n++){
				deltad += d[n]*d[n];
				
				alpha_num += - fPrimeVal[n]*d[n];
				alpha_denom += d[n]*fDoublePrimeVal[n]*d[n];
			}
			
			double alpha = alpha_num /alpha_denom;
			
			for(n=0; n<coordinateDim; n++){
				u[n] = u[n] + alpha*d[n];
			}
			j++;
			
			double alpha2_deltad = alpha*alpha*deltad;
			double epsilon2 = epsilon*epsilon;
			
			while(j<jmax && alpha2_deltad > epsilon2){
				double deltaOld=0;
				
				fPrimeVal = function.fPrime(u);
				for(n=0; n<coordinateDim; n++){
					r[n] = - fPrimeVal[n];
					deltaOld = deltaNew;
					
					deltaNew += r[n]*r[n];
				}
				
				double beta = deltaNew/ deltaOld;
				
				for(n=0; n<coordinateDim; n++){
					d[n] = r[n] + beta*d[n];
				}
				
				k++;
			}
			
			double rTd = 0;
			
			for(n=0; n<coordinateDim; n++){
				rTd += r[n]*d[n];
			}

			if(k==n || rTd <= 0){
				for(n=0; n<coordinateDim; n++){
					d[n] = r[n];
				}
				k=0;
			}
			function.getScalarEnergy();
			
			System.out.println("NonlinearCG within WHILE loop...");
			i++;
		}
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
	
}
