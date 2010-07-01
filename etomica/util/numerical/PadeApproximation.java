package etomica.util.numerical;

import Jama.Matrix;

/**
 * 
 * Class to approximate a power series with rational function (Pade Approximation)
 * 
 * 
 *  Pade[K/L] = A_K/ C_L = B_M
 *   with Kth, Lth and Mth order, where K + L = M  
 *  
 *  a0 + a1*x + a2*x^2 + ... + aK*x^K
 *  ----------------------------------  = b0 + b1*x + b2*x^2 + b3*x^3 + ... + b(K+L)*x^(K+L) [or bM*x^M]      
 *  c0 + c1*x + c2*x^2 + ... + cL*x^L
 * 
 *  - c0 is set to be 1
 * 
 * 
 * @author Tai Boon Tan
 *
 */
public class PadeApproximation {
	public PadeApproximation(double[] b, int K, int L){
		/*
		 * checking the order and we do not allow K+L > M
		 * K and L is the order of the power series of the numerator and denominator
		 *  or the rational function respectively
		 *  
		 * M is the order of the examined polynomial
		 */
		int M = (b.length -1);
		
		if((K+L) > M || (K-L)< 0.0){
			throw new RuntimeException("Pade Approximation class: K plus L should not exceed the M-order -OR- K-L should be greater than 0.0");
		}
		
		this.b = b;
		a = new double[K+1];
		c = new double[L+1];
		
		double[][] y = new double[L][L];
		double[] z = new double [L];
		
		for (int icol=0; icol<L; icol++){
			z[icol] = -b[K+1+icol];
		}
		
		for (int irow=0; irow<L; irow++){
			for (int icol=0; icol<L; icol++ ){
				if((K+irow-icol)< 0.0) {
					y[irow][icol] = 0.0;
				} else {
					y[irow][icol] = b[K+irow-icol]; 
				}
			}
		}
		
		Y = new Matrix(y);		
		Z = new Matrix(z, z.length);
		
	}
	
	public void solveCoefficients(){
		
		Matrix cMatrix = Y.solve(Z);
	
		/*
		 * determine "c" coefficients
		 */
		c[0] = 1.0;
		for (int i=0; i<c.length-1; i++){
			c[i+1] = cMatrix.get(i, 0);
		}
		
		/*
		 * determine "a" coefficients
		 */
		for (int i=0; i<a.length; i++){
			for (int j=0; j<c.length; j++){
				if((i-j)< 0.0 ) break;
				
				a[i] += b[i-j]*c[j];
				
			}
		}
	}
	
	public double[] getA(){
		return a;
	}
	
	public double[] getC(){
		return c;
	}
	
	public static void main(String[] args){
		
		double[] b = new double []{1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 2.0, 2.0};
		int K = 3;
		int L = 1;
		PadeApproximation pade = new PadeApproximation(b, K, L);
		pade.solveCoefficients();
		
		double[] aValues = pade.getA();
		double[] cValues = pade.getC();
		
		System.out.print("a:");
		for (int i=0; i<aValues.length; i++){
			System.out.print(" " + aValues[i]);
		}
		
		System.out.print("\nc:");
		for (int i=0; i<cValues.length; i++){
			System.out.print(" " + cValues[i]);
		}
		
	}
	
	protected double[] a, b, c;
	protected Matrix Y, Z;
}
