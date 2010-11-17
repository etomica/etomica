package etomica.util.numerical;

import etomica.util.FunctionMultiDimensionalDifferentiable;

public class ConjugateGradientMultiDimensional {
	
	/*
	 * 	Section 10.6 Numerical Derivative by Ridder's Method
	 * 	Numerical Recipes in FORTRAN, 2nd Edition, 1992
	 * 
	 * 	Given parameters:
	 * 
	 *  @ param  double[n] p  : n-dimensional points
	 *  		  function    : from FunctionMultiDimensionalDifferentiable which is n-Dimensional
	 *  			ftol      : convergence tolerance on the function value
	 *  
	 *  
	 * 	@ return   fret       : the value of function at the returned minimum location p
	 * 
	 * Author: Tai Boon Tan
	 * 					
	 */
	
	
	
	
	public ConjugateGradientMultiDimensional(){
		this.lineMinimization = new LineMinimization();
	}
	
	public void conjugateGradient(double[] p, double ftol, FunctionMultiDimensionalDifferentiable function){
		int j, its;
		int n = p.length;
		int ITMAX = 100;
		double EPS = 1E-10;
		
		this.minimumCoordinate = p;
		
		double gg, gam, fp, dgg;
		double[] g = new double[n];
		double[] h = new double[n];
		
		int[] derivative = new int[n];
		double[] df = new double[n];
		
		fp = function.f(p); 
		
		for (j=0; j<n; j++){
			
			derivative[j] = 1;
			df[j] = function.df(derivative, p);
			derivative[j] = 0;
		
		} 

		for (j=0; j<n; j++){
			g[j] = -df[j];
			df[j] = h[j] = g[j];
		}
		
		System.out.println("<ConjugateGradient>");
		for (its=1; its<ITMAX; its++){
			
			iteration = its;
			
			System.out.println("\n<CG> itereation# "+iteration);
			System.out.println("+++++++++++<ConjugateGradient> begin line minimization+++++++++");
			fret = lineMinimization.dLineMinimization(p, df, function);
			System.out.println("++++++++++<ConjugateGradient> after line minimization++++++++++++\n");
			minimumCoordinate = lineMinimization.getPcom();
//			System.out.println("x value is: "+ p[0]);
//			System.out.println("y value is: "+ p[1]);
			//System.out.println("Minimum value is: "+ fret);
			
			if (2.0*Math.abs(fret-fp) <= ftol*(Math.abs(fret)+Math.abs(fp)+EPS)){
				return;
			}
			
			fp = fret;
			
			for (j=0; j<n; j++){
				derivative[j] = 1;
				df[j] = function.df(derivative,p);
				derivative[j] = 0;
			}
			
			dgg = gg = 0.0;
			
			for (j=0; j<n; j++){
				gg += g[j]*g[j];
				dgg += (df[j]+g[j])*df[j];
				
			}
			
			if (gg == 0.0){
				//System.out.println("x value is: "+ p[0]);
				//System.out.println("y value is: "+ p[1]);
				//System.out.println("Minimum value is: "+ fret);
				return;
			}
			
			gam = dgg/gg;
			for (j=0; j<n; j++){
				g[j] = -df[j];
				df[j] = h[j] = g[j] + gam*h[j];
			}
			
		}
		
		throw new RuntimeException("Too many iteration in frprmn");
		
	}

	public double getFunctionMinimimumValue(){
		return this.fret;
	}
	
	public double[] getMinimumCoordinates(){
		return this.minimumCoordinate;
	}
	
	public int getNumIterations(){
		return this.iteration;
	}
	
	protected LineMinimization lineMinimization;
	protected double fret;
	protected int iteration;
	protected double[] minimumCoordinate;
	
}
