package etomica.conjugategradient;

import etomica.action.Activity;
import etomica.atom.AtomAgentManager;
import etomica.atom.IAtom;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.IVector;
import etomica.space.Space;

public class FiniteDifferenceDerivative {
	
	/*
	 * Section 5.7 Numerical Derivative by Ridder's Method
	 * Numerical Recipes in FORTRAN, 2nd Edition, 1992
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
	
	protected fFunction fFunction;
	protected IVector orientation ;
	protected double delta;
	protected double error;
	protected int dimension;
	
	
	public FiniteDifferenceDerivative(Box box, PotentialMaster potentialMaster){
		this.box = box;
		this.potentialMaster = potentialMaster;
		meterEnergy = new MeterPotentialEnergy(potentialMaster);
		allAtoms = new IteratorDirective();
		forceSum = new PotentialCalculationForceSum();
		
		MyAgentSource source = new MyAgentSource(box.getSpace());
		agentManager = new AtomAgentManager(source, box);
		forceSum.setAgentManager(agentManager);
	}
	
	public double[] finiteDerivative(fFunction function, double[] u, double h, int dimension, boolean hOptimizer){
		
		int coordinateDim = dimension;
		double[] dfridr = new double[coordinateDim]; 
		int ntab = 10;
		double con = 1.4;
		double con2 = con*con;
		double big = Math.pow(1, 30);
		double safe = 2.0;
		
		double errt, fac, hh;
		double[][][] a = new double[ntab][ntab][coordinateDim];
		double[] uPlus = new double[coordinateDim];
		double[] uMinus = new double[coordinateDim];
		
		hh = h;
		
		for (int p=0; p<coordinateDim; p++){  //loop over the p-th second derivatives
			
			for(int q=0; q<coordinateDim; q++){ // loop over the q-th generalized coordinate
				if(q==p){
					uPlus[q] = uPlus[q] + hh;
					uMinus[q] = uMinus[q] - hh;
				} else {
					uPlus[q] = u[q];
					uMinus[q] = u[q];
				}
			}
		
			a[0][0][p]= (function.fPrime(uPlus)[p] - function.fPrime(uMinus)[p])/(2.0*hh);
		
			if (!hOptimizer) {
				dfridr[p] = a[0][0][p];
				continue;
			}
			
			double err = big;
			
			for(int i=1; i<ntab; i++){
				hh = hh /con;
				a[0][i][p] = (function.fPrime(uPlus)[p] - function.fPrime(uMinus)[p])/(2.0*hh);
				fac = con2;
				
				for(int j=1; j<i; j++){
					a[j][i][p] = (a[j-1][i][p]*fac - a[j-1][i-1][p])/(fac-1);
					fac = con2*fac;
					errt = Math.max(Math.abs(a[j][i][p]-a[j-1][i][p]), Math.abs(a[j][i][p]-a[j-1][i-1][p]));
					
					if (errt <= err){
						err = errt;
						dfridr[p] = a[j][i][p];
					}
				}
				
				if (Math.abs(a[i][i][0]-a[i-1][i-1][0]) >= safe*err){
					break;
				}
			}
		
		} //end of looping p
		
		return dfridr;
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
