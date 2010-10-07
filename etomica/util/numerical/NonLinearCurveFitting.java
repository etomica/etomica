package etomica.util.numerical;

import Jama.LUDecomposition;
import Jama.Matrix;

/**
 * 
 * 
 * @author taitan
 *
 */
public class NonLinearCurveFitting {
	public NonLinearCurveFitting(String filename){

		funcs = new FittingFunctionNonLinear();
		
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
		
		ma = 2;
		atry = new double[ma];
		da = new double[ma];
		
		covar = new double[ma][ma];
		a = new double[]{1.0, 1.0};
		beta = new double[ma];
		alpha = new double[ma][ma];
		alamda = -1.0;
		
	}
	
	public double[] findParameter(){

		int counter =0;
		while(iterate){
			mrqmin();
			System.out.println("**"+a[0]+" "+a[1]+" " + chisq);
			++counter;
			if(counter == 100){
				System.out.println("SCREW!");
				iterate = false;
			}
		}
		
		alamda=0.0;
		//mrqmin();
		
		return a;
	}
	
	public void mrqmin(){
	
		if(alamda < 0.0){
			alamda = 0.01;
		
			mrqcofA(a);
			ochisq = chisq;
			for (int j=0; j<ma; j++){
				atry[j] = a[j];
			}
		}
		
//			System.out.println("\na: "+a[0]+" "+a[1]);
//			System.out.println("chi: "+chisq+" "+ochisq);
//			if(Double.isInfinite(alamda)){
//				System.exit(1);
//			}
//			System.out.println("alamda: " + alamda);
			
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
//			System.out.println("@@ aMatrix");
//			for (int irow=0; irow<aMatrix.getRowDimension(); irow++){
//				for (int icol=0; icol<aMatrix.getColumnDimension(); icol++){
//					System.out.print(aMatrix.get(irow, icol)+" ");
//				}	
//				System.out.println();
//			}
//			System.out.println("@@ bMatrix");
//			for (int irow=0; irow<bMatrix.getRowDimension(); irow++){
//				for (int icol=0; icol<bMatrix.getColumnDimension(); icol++){
//					System.out.print(bMatrix.get(irow, icol)+" ");
//				}	
//				System.out.println();
//			}
			
			LUDecomposition lu = new LUDecomposition(aMatrix);
			Matrix xMatrix = lu.solve(bMatrix);
			da = xMatrix.getColumnPackedCopy();
//			System.out.println("da: " + da[0] + " " + da[1]);
			
//			if(alamda == 0.0){
//				
//				covsrt(mfit, covar);
//				return;
//			}
			
			System.out.println("a: " + a[0] + " " + a[1]);
			for (int l=0; l<ma; l++){
				atry[l] = a[l] + da[l];
			}
			mrqcof(atry);
			
			if(chisq < ochisq){
				System.out.println("**** newChiSq < oldChiSq ****  " + chisq + " " + ochisq );
		
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
			
				//System.exit(1);
			} else {
				System.out.println("++++++ newChiSq > oldChiSq ++++++  " + chisq + " " + ochisq );
				alamda *= 10;
				chisq = ochisq;
			
			}
			return;
	}
	
	
	/*
	 * calculate alpha[][], beta[] and chisq 
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
	
	public void mrqcofA(double[] atry){
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
			double dy = y[i] - funcs.f(atry, x[i]); 
			
			for (int j=0; j<ma; j++){
				dyda[j] = funcs.df(j, atry, x[i]);
				wt = dyda[j]*sig2i;
				
				for (int k=j; k<ma; k++){
					alpha[j][k] += wt*funcs.df(k, atry, x[i]);
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
		String filename = "/tmp/foo";
		
		NonLinearCurveFitting nonLinFit = new NonLinearCurveFitting(filename);
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
