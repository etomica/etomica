package etomica.conjugategradient;

import etomica.util.FunctionMultiDimensional;
import etomica.util.FunctionMultiDimensionalDifferentiable;

public class FiniteDifferenceDerivativeWOBox implements FunctionMultiDimensionalDifferentiable{
	
	/*
	 * Section 5.7 Numerical Derivative by Ridder's Method
	 * Numerical Recipes in FORTRAN, 2nd Edition, 1992
	 * 
	 * @author Tai Tan
	 */
	
	
	protected FunctionMultiDimensional fFunction;
	protected double error;
	protected double h;
	protected boolean hOptimizer;
	
	public FiniteDifferenceDerivativeWOBox(FunctionMultiDimensional fFunction){

		this.fFunction = fFunction;
		h = 0.00001;
		hOptimizer = false;
	}
	
	public double function(double[] u){
		return fFunction.function(u);
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
		
			a[0][0][p]= (fFunction.function(uPlus) - fFunction.function(uMinus))/(2.0*hh);
			System.out.println(" a[0][0]["+p+"] is: "+a[0][0][p]);
			if (!hOptimizer) {
				dfdx[p] = a[0][0][p];
				continue;
			}
			
			double err = big;
			
			for(int i=1; i<ntab; i++){
				hh = hh /con;
				a[0][i][p] = (fFunction.function(uPlus) - fFunction.function(uMinus))/(2.0*hh);
				System.out.println(" a[0]["+i+"]["+p+"] is: "+a[0][i][p]);
				fac = con2;
				
				for(int j=1; j<i; j++){
					a[j][i][p] = (a[j-1][i][p]*fac - a[j-1][i-1][p])/(fac-1);
					System.out.println(" a["+j+"]["+i+"]["+p+"] is: "+a[j][i][p]);
					fac = con2*fac;
					errt = Math.max(Math.abs(a[j][i][p]-a[j-1][i][p]), Math.abs(a[j][i][p]-a[j-1][i-1][p]));
					//System.out.println("errt is: "+errt);
					
					if (errt <= err){
						err = errt;
						dfdx[p] = a[j][i][p];
						System.out.println("in errt<= err, dfdx["+p+"] is: "+a[j][i][p]);
					}
				}
				
				if (Math.abs(a[i][i][p]-a[i-1][i-1][p]) >= safe*err){
					break;
				}
			}
		
		} //end of looping p
		
		return dfdx;
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
