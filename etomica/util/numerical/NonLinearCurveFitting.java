package etomica.util.numerical;

import Jama.LUDecomposition;
import Jama.Matrix;

/**
 * 
 * Class that compute 
 * 
 * Reference: NUMERICAL RECIPES in Fortran 2nd Ed.
 * 			  Levenberg-Marquardt Method (Chapter 15)
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
 * - Return parameters, a[]
 * 
 * Note: Eqs. are from Reference
 * 
 * @author Tai Boon Tan
 *
 */
public class NonLinearCurveFitting {
	public NonLinearCurveFitting(String filename, int M, int N){

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
		
		this.ma = M+1+N;
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
		funcs = new FittingFunctionNonLinear(M, N);
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
		
		// To determine the covariance matrix once the parameters
		// have converged.
		if(alamda == 0.0){
			return;
		}
			
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
	
	public double computeRMSD(double[] a){
		double RMSD = 0;
		double diff;
		
		for(int i=0; i<ndata; i++){
			diff = ((funcs.f(a, x[i])-y[i])/sig[i]);
			RMSD += diff*diff;
			
		}
		
		return Math.sqrt(RMSD/ndata);
	}
	
	public static void main(String[] args){
		String filename = "/tmp/foo_stat";
		int M = 0;
		int N = 1;
		if(args.length > 0){
			filename = args[0];
		}
		if(args.length > 1){
			M = Integer.parseInt(args[1]);
		}
		if(args.length > 2){
			N = Integer.parseInt(args[2]);
		}
		
		
		NonLinearCurveFitting nonLinFit = new NonLinearCurveFitting(filename, M, N);
		double[] a = nonLinFit.findParameter();
		for(int i=0; i<M; i++){
			System.out.println("a[" + i +"]: "+a[i]);
		}
		System.out.println("K: " + a[M]);
		for(int i=M+1; i<M+1+N; i++){
			System.out.println("b[" + (i-(M+1)) +"]: "+a[i]);
		}
		
		double RMSD = nonLinFit.computeRMSD(a);
		System.out.println("RMSD: " + RMSD);
		
		//Covariance Matrix
//		for(int i=0; i<nonLinFit.covar.length; i++){
//			for(int j=0; j<nonLinFit.covar[0].length; j++){
//				System.out.print(nonLinFit.covar[i][j]+" ");
//			}	
//			System.out.println();
//		}
	}
	
	protected FittingFunctionNonLinear funcs;
	protected int ma, ndata;
	protected double[] a, beta, sig, x, y;
	protected double[][] covar, alpha;
	protected boolean iterate = true;
	protected double chisq, ochisq, alamda;
	protected double[] atry, da;
}
