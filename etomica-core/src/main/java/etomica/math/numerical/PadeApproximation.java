/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.numerical;

import Jama.Matrix;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

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
	public PadeApproximation(double[] c, int N, int D){
		/*
		 * checking the order and we do not allow K+L > M
		 * K and L is the order of the power series of the numerator and denominator
		 *  or the rational function respectively
		 *  
		 * M is the order of the examined polynomial
		 */
		int M = (c.length -1);
		
		if((N+D) > M || D<0){
		    System.err.println("N+D="+(N+D)+"  and  M="+M);
			throw new RuntimeException("Pade Approximation class: N plus D should not exceed the M-order -OR- D < 0");
		}
		
		this.c = c;
		a = new double[N+1];
		b = new double[D+1];
		
		double[][] y = new double[D][D];
		double[] z = new double [D];
		
		for (int icol=0; icol<D; icol++){
			z[icol] = -c[N+1+icol];
		}
		
		for (int irow=0; irow<D; irow++){
			for (int icol=0; icol<D; icol++ ){
				if((N+irow-icol)< 0.0) {
					y[irow][icol] = 0.0;
				} else {
					y[irow][icol] = c[N+irow-icol]; 
				}
			}
		}
		
		Y = D==0 ? new Matrix(0,0) : new Matrix(y);		
		Z = new Matrix(z, z.length);
		
	}
	
	public void solveCoefficients(){
		
		Matrix bMatrix = Y.solve(Z);
	
		/*
		 * determine "c" coefficients
		 */
		b[0] = 1.0;
		for (int i=0; i<b.length-1; i++){
			b[i+1] = bMatrix.get(i, 0);
		}
		
		/*
		 * determine "a" coefficients
		 */
		for (int i=0; i<a.length; i++){
			for (int j=0; j<b.length; j++){
				if((i-j)< 0.0 ) break;
				
				a[i] += c[i-j]*b[j];
				
			}
		}
	}
	
	public double[] getA(){
		return a;
	}
	
	public double[] getB(){
		return b;
	}
	
	public static void main(String[] args){
		
		VirialParam params = new VirialParam();
	    ParseArgs.doParseArgs(params, args);
		
		double[] c = params.c;
				
		int n = params.n;
		int d = params.d;
		if (n+d != c.length-1) {
		    throw new RuntimeException("n+d must equal the number of coefficients");
		}
		PadeApproximation pade = new PadeApproximation(c, n, d);
		pade.solveCoefficients();
		
		double[] aValues = pade.getA();
		double[] bValues = pade.getB();
		
		for (int i=0; i<aValues.length; i++){
			System.out.println("a " + i +"  "+ aValues[i]);
		}
		
		System.out.println();
		for (int i=0; i<bValues.length; i++){
			System.out.println("b " + i +"  " + bValues[i]);
		}
	}
	
	protected double[] a, b, c;
	protected Matrix Y, Z;
	
	public static class VirialParam extends ParameterBase {
	    public double[] c = new double[]{     1.0, 3.712218666, 5.55200, 
                    1.44261,     -1.6883,  1.8935, 
                    -1.700,        0.44,  3.0589};
	    public int d = 4, n = 4;
	}
}
