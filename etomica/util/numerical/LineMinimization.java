package etomica.util.numerical;

import etomica.util.FunctionMultiDimensionalDifferentiable;

public class LineMinimization {
	
	/*
	 * 	Section 10.6 Numerical Derivative by Ridder's Method
	 * 	Numerical Recipes in FORTRAN, 2nd Edition, 1992
	 * 
	 * 	Given parameters:
	 * 
	 *  @ param  double[n] p  : n-dimensional points
	 *  		 double[n] xi : n-dimensional directions			
	 *  		  function    : from FunctionMultiDimensionalDifferentiable which is n-Dimensional
	 *  	
	 *  moves and resets p to where the function function(p) takes on a minimum along the direction xi from p,
	 *  and replaces xi by the actual vector displacement that p was moved.	And:
	 *  
	 * 	@ return   function(p)  : the value of function at the returned location p, which is dbrent[1]
	 * 
	 */
	
	
	
	public LineMinimization(){
		this.f1dim = new Function1d();
		this.bisectionMethodMinimizationBracket = new BisectionMethodMinimumBracket();
		this.brentMethodwDerivative = new BrentMethodwDerivative();
	}
	

	public double dLineMinimization(double[] p, double[] xi, FunctionMultiDimensionalDifferentiable function){
		double TOL = 0.0002;
		
		int j;
		double xx, xmin, bx, ax;
		int n = p.length;
		ncom = n;
		pcom = new double[n];
		xicom = new double[n];
		
		
		for(j=0; j<n; j++){
			pcom[j] = p[j];
			xicom[j] = xi[j];
		}
		
		ax = 0.0;
		xx = 1.0;
		bx = Double.NaN;
		
		f1dim.setFunction(function);
		f1dim.setPcom(pcom);
		f1dim.setXicom(xicom);

		double[] x = new double[] {ax, xx, bx}; 
		bisectionMethodMinimizationBracket.mnbrak(x, f1dim);
		
		double[] dbrentValue = brentMethodwDerivative.dbrent(x, f1dim, TOL);
		xmin = dbrentValue[0];
		
		for (j=0; j<n; j++){
			xi[j] *= xmin;
			p[j] += xi[j];
		}
	
		return dbrentValue[1];
	}
	
	protected int ncom;
	protected double[] pcom, xicom, nrfunc;
	protected BisectionMethodMinimumBracket bisectionMethodMinimizationBracket;
	protected BrentMethodwDerivative brentMethodwDerivative;
	protected LineMinimization lineMinimization;
	public Function1d f1dim;
	
	
	
}
