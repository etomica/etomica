package etomica.util.numerical;

import Jama.LUDecomposition;
import Jama.Matrix;

/**
 * 
 * Class that compute 
 * 
 * Reference: NUMERICAL RECIPES in Fortran 2nd Ed.
 * 			  Levenberg-Mrquardt Method (Chapter 15)
 * 	          - combination of the use of 
 * 				(a) Inverse-Hessian Method  [Eq. 15.5.9]
 *              (b) Steepest Descent Method [Eq. 15.5.10]
 *              
 * Adapting codes from subroutines: mrqmin, covsrt, and mrqcof
 * 
 * The program:
 * - Takes in user-defined function FittingFunctionNonLinear
 *   where the func and derivative are defined
 * - Calculates for beta[] and alpha[][] as defined in [Eq 15.5.8 through 15.5.11]
 * - Solves for da [Eq 15.5.14] 
 * - Return paramters, a[]
 * 
 * Note: Eqs. are from Reference
 * 
 * @author Tai Boon Tan
 *
 */
public class NonLinearCurveFitting {
	public NonLinearCurveFitting(String filename, int ma){

		// Reading data values from file [x, y, sigma]
		double[][] value = ArrayReader1D.getFromFile(filename);
		ndata =value.length;
		x = new double[ndata];
		y = new double[ndata];
		sig = new double[ndata];

		for (int i=0; i<ndata; i++){
			x[i]   = value[i][0];
			y[i]   = value[i][1];
			sig[i] = value[i][2];
		}
		
		this.ma = ma;
		a = new double[ma];
		atry = new double[ma];
	
		da = new double[ma];	
		beta = new double[ma];
		
		covar = new double[ma][ma];
		alpha = new double[ma][ma];
		
		// Give initial guesses for the parameter
		for (int i=0; i<a.length; i++){
			a[i] = 1.0;
		}

		alamda = -1.0;
		funcs = new FittingFunctionNonLinear();
	}
	
	public double[] findParameter(){

		while(iterate){
			mrqmin();
		}
		
		alamda=0.0;
		mrqmin();
		
		return a;
	}
	
	/*
	 * Marquardt's Method (one iteration)
	 */
	public void mrqmin(){
	
		if(alamda < 0.0){
			alamda = 0.01;
		
			mrqcofA(a);
			ochisq = chisq;
			for (int j=0; j<ma; j++){
				atry[j] = a[j];
			}
		}
			
		for(int j=0; j<ma; j++){
			for(int k=0; k<ma; k++){
				if(j==k){
					covar[j][k] = alpha[j][k]*(1.0+alamda);	
				} else {
					covar[j][k] = alpha[j][k];
				}
			}
			da[j] = beta[j];
		}
			
		/*
		 * solving deltaA; Equation (15.5.14)
		 */
		Matrix aMatrix = new Matrix(covar);
		Matrix bMatrix = new Matrix(da, da.length);
			
		LUDecomposition lu = new LUDecomposition(aMatrix);
		Matrix xMatrix = lu.solve(bMatrix);
		da = xMatrix.getColumnPackedCopy();
		
//		if(alamda == 0.0){
//			covsrt(mfit, covar);
//			return;
//		}
			
		for (int l=0; l<ma; l++){
			atry[l] = a[l] + da[l];
		}
	
		mrqcof(atry);
		
		if(chisq < ochisq){
	
			if((Math.abs(chisq-ochisq)<1e-2)){
				iterate = false;
			}
			alamda *= 0.1;
			ochisq = chisq;
			
			for(int j=0; j<ma; j++){
				for (int k=0; k<ma; k++){
					alpha[j][k] = covar[j][k];
				}
				beta[j] = da[j];
				a[j] = atry[j];
			}
			
		} else {
			alamda *= 10;
			chisq = ochisq;
			
		}
	}
	
	
	/*
	 * Calculate alpha[][], beta[] and chisq 
	 */
	public void mrqcof(double[] atry){
		double wt;
		double [] dyda = new double[ma];
		
		// Initialize (symmetric) alpha, beta
		for (int j=0; j<ma; j++){
			for (int k=0; k<ma; k++){
				covar[j][k] = 0.0;
			}
			da[j] = 0.0;
		}
		chisq = 0.0;
		
		// Summation loop over all data
		for(int i=0; i<ndata; i++){
			double sig2i = 1.0/(sig[i]*sig[i]);
			double dy = y[i] - funcs.f(atry, x[i]); 
			
			for (int j=0; j<ma; j++){
				dyda[j] = funcs.df(j, atry, x[i]);
				wt = dyda[j]*sig2i;
				
				for (int k=j; k<ma; k++){
					covar[j][k] += wt*funcs.df(k, atry, x[i]);
				}
				da[j] += dy*wt;
			}
			// find chi^2
			chisq += dy*dy*sig2i;
		}
		
		// Fill in the symmetric side
		for(int j=0; j<ma; j++){
			for (int k=j+1; k<ma; k++){
				covar[k][j] = covar[j][k];
			}
		}
	}
	
	public void mrqcofA(double[] a){
		double wt;
		double [] dyda = new double[ma];
		
		// Initialize (symmetric) alpha, beta
		for (int j=0; j<ma; j++){
			for (int k=0; k<ma; k++){
				alpha[j][k] = 0.0;
			}
			beta[j] = 0.0;
		}
		chisq = 0.0;
		
		// Summation loop over all data
		for(int i=0; i<ndata; i++){
			double sig2i = 1.0/(sig[i]*sig[i]);
			double dy = y[i] - funcs.f(a, x[i]); 
			
			for (int j=0; j<ma; j++){
				dyda[j] = funcs.df(j, a, x[i]);
				wt = dyda[j]*sig2i;
				
				for (int k=j; k<ma; k++){
					alpha[j][k] += wt*funcs.df(k, a, x[i]);
				}
				beta[j] += dy*wt;
			}
			// find chi^2
			chisq += dy*dy*sig2i;
		}
		
		// Fill in the symmetric side
		for(int j=0; j<ma; j++){
			for (int k=j+1; k<ma; k++){
				alpha[k][j] = alpha[j][k];
			}
		}
	}
//	public void covsrt(int mfit, double[][] covar){
//		for (int i=mfit+1; i<ma; i++){
//			for (int j=0; j<i; j++){
//				covar[i][j] = 0.0;
//				covar[j][i] = 0.0;
//			}
//		}
//		
//		double swap;
//		int k=mfit-1;
//		
//		for (int j=ma-1; j>=0; j--){
//			if(ia[j] != 0){
//				for (int i=0; i<ma; i++){
//					swap = covar[i][k];
//					covar[i][k] = covar[i][j];
//					covar[i][j] = swap;
//					
//				}
//				
//				for (int i=0; i<ma; i++){
//					swap = covar[k][i];
//					covar[k][i] = covar[j][i];
//					covar[j][i] = swap;
//					
//				}
//				k -= 1;
//			}
//		}
//	}
	
	public static void main(String[] args){
		String filename = "/tmp/foo_stat";
		int ma = 2;
		NonLinearCurveFitting nonLinFit = new NonLinearCurveFitting(filename, ma);
		double[] a = nonLinFit.findParameter();
		for(int i=0; i<a.length; i++){
			System.out.println("a[" + i +"]: "+a[i]);
		}
	}
	
	protected FittingFunctionNonLinear funcs;
	protected int ma, ndata;
	protected double[] a, beta, sig, x, y;
	protected double[][] covar, alpha;
	protected boolean iterate = true;
	protected double chisq, ochisq, alamda;
	protected double[] atry, da;
}
