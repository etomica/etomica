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
	 */
	
	
	
	
	public ConjugateGradientMultiDimensional(){
		this.lineMinimization = new LineMinimization();
	}
	
	public void conjugateGradient(double[] p, double ftol, FunctionMultiDimensionalDifferentiable function){
		int j, its;
		int n = p.length;
		int ITMAX = 300;
		double EPS = Math.pow(10, -10);
		
		this.minimumCoordinate = p;
		
		double gg, gam, fp, dgg;
		double[] g = new double[n];
		double[] h = new double[n];
		
		int[] derivative = new int[n];
		double[] df = new double[n];
		
		fp = function.f(p); 
		
		for (j=0; j<n; j++){
			
			derivative[j] = 1;
			
			df[j] = function.df(derivative,p);
			System.out.println("The d["+j+"] is: "+ df[j]);
			derivative[j] = 0;
		}
		
		for (j=0; j<n; j++){
			g[j] = -df[j];
			df[j] = h[j] = g[j];
		}
		
		for (its=0; its<ITMAX; its++){
			
			this.iteration = its;
			
			for (int i=0; i<p.length; i++){
				System.out.println("The number of iteration: "+its);
				System.out.println("u["+i+"] is: "+p[i]);
			}
			
			fret = lineMinimization.dLineMinimization(p, df, function);
			
			//System.out.println("x value is: "+ p[0]);
			//System.out.println("y value is: "+ p[1]);
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
	
	public double getNumIterations(){
		return this.iteration;
	}
	
	protected LineMinimization lineMinimization;
	protected double fret;
	protected int iteration;
	protected double[] minimumCoordinate;
	
}
