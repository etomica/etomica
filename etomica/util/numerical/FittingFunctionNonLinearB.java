/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.util.numerical;


/**
 * 
 * This class is coded for specific nonlinear model that decribes the Equation:
 *  f(x) = [ 1 + exp(-Kx^L) ] * polynomialA
 * 
 * M is the input parameter that determine the order of polynomialA
 *  
 * polynomialA is given as: 
 *                 M
 *  polynomialA = sum [ a_m * x^m ]
 *                m=1
 * 
 * for example, when M=3; polynomialA = a_1 * x + a_2 * x^2 + a_3 * x^3
 * 
 * 
 * @author Tai Boon Tan
 *
 */
public class FittingFunctionNonLinearB{
	public FittingFunctionNonLinearB(int M){
		this.M = M;
	}

	public double f(double[] a, double x) {
		double sumPolyA = 0.0;
		double expKxL = Math.exp(a[0]*Math.pow(x,a[1]));

		//polynomialA
		double xMultiply = 1;
		for(int i=2; i<M+2; i++){
			xMultiply *= x;
			sumPolyA += a[i]*xMultiply;
			
		}
		
		return (1 + expKxL )* sumPolyA;
		
	}

	/*
	 * d is derivatives
	 * a is array of parameters
	 * x is the x-value
	 */
	public double df(int d, double[] a, double x) {
		double dsumPolyA = 0.0;
		double expKxL = Math.exp(a[0]*Math.pow(x, a[1]));
		

		double xMultiply = 1;
		if(d >= 2 ){
			for (int i=2; i<a.length; i++){
				xMultiply *= x;
				if(d==i){
					return (1 + expKxL) * xMultiply;
				}
			}
		} else {
		
			//polynomialA
			xMultiply = 1;
			for(int i=2; i<M+2; i++){
				xMultiply *= x;
				dsumPolyA += a[i]*xMultiply;
				
			}
			
			if(d==0){
				return Math.pow(x, a[1])*expKxL*dsumPolyA;
				
			} else if (d==1){
				return a[0]*Math.log(x)*Math.pow(x, a[1])*expKxL * dsumPolyA;
				
			}
		
		}
		System.out.println("screw!");
		return Double.NaN;

	}
	
	protected int M;
}
