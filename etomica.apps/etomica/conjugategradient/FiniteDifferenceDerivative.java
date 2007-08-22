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
import etomica.space.Space;
import etomica.util.FunctionMultiDimensional;
import etomica.util.FunctionMultiDimensionalDifferentiable;

public class FiniteDifferenceDerivative implements FunctionMultiDimensionalDifferentiable{
	
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
	
	protected FunctionMultiDimensional fFunction;
	protected double error;
	protected double h;
	protected boolean hOptimizer;
	
	public FiniteDifferenceDerivative(Box box, PotentialMaster potentialMaster, FunctionMultiDimensional fFunction){
		this.box = box;
		this.potentialMaster = potentialMaster;
		this.fFunction = fFunction;
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
		return fFunction.f(u);
	}
	
	public double[] dfdx(double[] u){
		
		int coordinateDim = u.length;
		double[] dfdx = new double[coordinateDim]; 
		int ntab = 10;
		double con = 1.4;
		double con2 = con*con;
		double big = Math.pow(10, 30);
		double safe = 2.0;
		
		double errt, fac, hh;
		double[][][] a = new double[ntab][ntab][coordinateDim];
		double[] uPlus = new double[coordinateDim];
		double[] uMinus = new double[coordinateDim];
		
		hh = h;
		
		for (int p=0; p<coordinateDim; p++){  //loop over the p-times second derivatives
			
			for(int q=0; q<coordinateDim; q++){ // loop over the q-times generalized coordinate
				if(q==p){
					uPlus[q] = u[q] + hh;
					uMinus[q] = u[q] - hh;
				} else {
					uPlus[q] = u[q];
					uMinus[q] = u[q];
				}
			}
		
			a[0][0][p]= (fFunction.f(uPlus) - fFunction.f(uMinus))/(2.0*hh);
		
			if (!hOptimizer) {
				dfdx[p] = a[0][0][p];
				continue;
			}
			
			double err = big;
			
			for(int i=1; i<ntab; i++){
				hh = hh /con;
				a[0][i][p] = (fFunction.f(uPlus) - fFunction.f(uMinus))/(2.0*hh);
				fac = con2;
				
				for(int j=1; j<i; j++){
					a[j][i][p] = (a[j-1][i][p]*fac - a[j-1][i-1][p])/(fac-1);
					fac = con2*fac;
					errt = Math.max(Math.abs(a[j][i][p]-a[j-1][i][p]), Math.abs(a[j][i][p]-a[j-1][i-1][p]));
					
					if (errt <= err){
						err = errt;
						dfdx[p] = a[j][i][p];
					}
				}
				
				if (Math.abs(a[i][i][p]-a[i-1][i-1][p]) >= safe*err){
					break;
				}
			}
		
		} //end of looping p
		
		return dfdx;
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
