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
	
	public double finiteDerivative(fFunction function, IVector x, double h){
		
		double dfridr=1; //for the moment
		int ntab = 10;
		double con = 1.4;
		double con2 = con*con;
		double big = Math.pow(1, 30);
		double safe = 2.0;
		
		double errt, fac, hh;
		double[][] a = new double[ntab][ntab];
		
		
		hh = h;
		// a[0][0] = (function(x+hh)-function(x-hh))/(2.0*hh)
		double err = big;
		
		for(int i=1; i<ntab; i++){
			hh = hh /con;
			// a[0][i] = (function(x+hh)-function(x-hh))/(2.0*hh)
			fac = con2;
			
			for(int j=1; j<i; j++){
				a[j][i] = (a[j-1][i]*fac - a[j-1][i-1])/(fac-1);
				fac = con2*fac;
				errt = Math.max(Math.abs(a[j][i]-a[j-1][i]), Math.abs(a[j][i]-a[j-1][i-1]));
				
				if (errt <= err){
					err = errt;
					dfridr = a[j][i];
				}
			}
			
			if (Math.abs(a[i][i]-a[i-1][i-1]) >= safe*err){
				return dfridr;
			}
		}
		
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
