/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.util.numerical;

import Jama.LUDecomposition;
import Jama.Matrix;

/**
 * Similar as <NonLinearCurveFitting> but fit the difference
 * 
 * See <NonLinearCurveFitting> for fundamental implementation details
 * 
 * @author Tai Boon Tan
 *
 */
public class NonLinearCurveDiffFitting {
	public NonLinearCurveDiffFitting(String filename, int M, int N){

		// Reading data values from file [x, y, sigma]
		double[][] value = ArrayReader1D.getFromFile(filename);
		ndata =value.length;
		x0 = new double[ndata];
		x1 = new double[ndata];
		y = new double[ndata];
		sig = new double[ndata];

		for (int i=0; i<ndata; i++){
			x0[i]  = value[i][0];
			x1[i]  = value[i][1];
			y[i]   = value[i][2];
			sig[i] = value[i][3];
		}
		
		this.ma = M+2+N;
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
			//System.out.println("chisq: " + chisq);
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
			double dy = y[i] - (funcs.f(atry, x1[i]) - funcs.f(atry, x0[i]) ); 
			
			for (int j=0; j<ma; j++){
				dyda[j] = (funcs.df(j, atry, x1[i]) - funcs.df(j, atry, x0[i]));
				wt = dyda[j]*sig2i;
				
				for (int k=j; k<ma; k++){
					covar[j][k] += wt*  (funcs.df(k, atry, x1[i]) - funcs.df(k, atry, x0[i]));
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
			double dy = y[i] - (funcs.f(a, x1[i]) - funcs.f(a, x0[i])); 
			
			for (int j=0; j<ma; j++){
				dyda[j] = (funcs.df(j, a, x1[i]) - funcs.df(j, a, x0[i]));
				wt = dyda[j]*sig2i;
				
				for (int k=j; k<ma; k++){
					alpha[j][k] += wt*(funcs.df(k, a, x1[i]) - funcs.df(k, a, x0[i]));
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
			diff = ((funcs.f(a, x1[i])  -funcs.f(a, x0[i])) -y[i])/sig[i];
			RMSD += diff*diff;
			
		}
		
		return Math.sqrt(RMSD/ndata);
	}
	
	public static void main(String[] args){
		String filename = "Adiff10.dat";
		int M = 3;
		int N = 10;
		if(args.length > 0){
			filename = args[0];
		}
		if(args.length > 1){
			M = Integer.parseInt(args[1]);
		}
		if(args.length > 2){
			N = Integer.parseInt(args[2]);
		}
		
		
		NonLinearCurveDiffFitting nonLinFit = new NonLinearCurveDiffFitting(filename, M, N);
		double[] a = nonLinFit.findParameter();
//		for(int i=0; i<M; i++){
//			System.out.println("a[" + (i+1) +"]: "+a[i]);
//		}
//		System.out.println("K: " + a[M]);
//		System.out.println("L: " + a[M+1]);
//		for(int i=M+2; i<M+2+N; i++){
//			System.out.println("b[" + (i-(M+1)) +"]: "+a[i]);
//		}
//		
//		double RMSD = nonLinFit.computeRMSD(a);
//		System.out.println("RMSD: " + RMSD);
		
		for (int i=0; i<a.length; i++){
			System.out.printf("%3.10e", a[i]);
			System.out.print(" ");
		}
		
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
	protected double[] a, beta, sig, x0, x1, y;
	protected double[][] covar, alpha;
	protected boolean iterate = true;
	protected double chisq, ochisq, alamda;
	protected double[] atry, da;
}
