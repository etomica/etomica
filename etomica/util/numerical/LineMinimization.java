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
		double TOL = 2.1e-9;
		/*
		 * xi is the first derivative
		 */
		
		double xmin;
		int n = p.length;
		ncom = n;
		pcom = new double[n];
		xicom = new double[n];

		for(int j=0; j<n; j++){
			pcom[j] = p[j];
			xicom[j] = xi[j];
			System.out.println("<lineMinimize>: pcom: " + pcom[j] + " ;xicom: " + xicom[j]);
		}
		
		// Values for bracketing
		double ax = -1.5e-8;
		double xx =  1.5e-8;
		double bx = Double.NaN;
		
		f1dim.setFunction(function);
		f1dim.setPcom(pcom);
		f1dim.setXicom(xicom);

		double[] x = new double[] {ax, xx, bx}; 
		
		System.out.println("\n<lineMinimize> BEGIN ******BRACKETING****");
		x=bisectionMethodMinimizationBracket.mnbrak(x, f1dim);
		System.out.println("<lineMinimize> END ******BRACKETING****");
		System.out.println("bracketing the minimum: " +x[0]+" "+x[1]+" "+x[2]);
//		System.out.println("exit in LineMinimization");
//		System.exit(1);
		System.out.println("\n<lineMinimize> BEGIN *******dbrent******");
		double[] dbrentValue = brentMethodwDerivative.dbrent(x, f1dim, TOL);
		System.out.println("<lineMinimize> END *******dbrent******");
		System.out.println("dbrentvalue: " + dbrentValue[0] +" "+dbrentValue[1]);
	

		xmin = dbrentValue[0];
		System.out.println("\n<lineMinimize>*******for loop*******");
		System.out.println("<lineMinimize> xmin: " + xmin);
		for (int j=0; j<n; j++){
			xi[j] *= xmin;
			p[j] += xi[j];
			//System.out.println("<lineMinimize>#########LineMinimization: xi: " + xi[j] + " ;p[j]: " + p[j]);
		}
	
		for (int j=0; j<n; j++){
			System.out.print(p[j] + ", ");
			if (j>0&&(j+1)%5==0){
				System.out.println("");
			}
		}
		System.out.println("Lowest Energy found: "+ dbrentValue[1]);
		return dbrentValue[1];
	}
	
	protected int ncom;
	protected double[] pcom, xicom, nrfunc;
	protected BisectionMethodMinimumBracket bisectionMethodMinimizationBracket;
	protected BrentMethodwDerivative brentMethodwDerivative;
	protected LineMinimization lineMinimization;
	public Function1d f1dim;
	
	
	
}
