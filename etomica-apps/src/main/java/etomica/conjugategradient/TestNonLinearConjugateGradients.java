/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.conjugategradient;

import etomica.math.function.FunctionMultiDimensionalDifferentiable;
import etomica.math.numerical.FiniteDifferenceDerivative;

public class TestNonLinearConjugateGradients implements FunctionMultiDimensionalDifferentiable{
	
	public TestNonLinearConjugateGradients(){
		
	}
	
	public double f(double[] u){
		return u[0]*u[0] + u[1]*u[1];
	}
	
	public double df(int[] d, double[] u){
		
		if (d[0]==0){
			if (d[1]==0){
				return f(u);
			} else if (d[1]==1){
				return 2*u[1];
			} else if (d[1]==2){
				return 2;
			} 
		} else if (d[0]==1){
			return(d[1]==0)? 2*u[0]: 0.0;
		} else if (d[0]==2){
			return(d[1]==0)? 2: 0.0;
		}
		return 0;
		
	}
	
	public int getDimension(){
		return 1;
	}
	
	public static void main(String args[]){
		double[] uNew = new double[] {.2, .2};
		TestNonLinearConjugateGradients testFunction = new TestNonLinearConjugateGradients();
		
		FiniteDifferenceDerivative finiteDifferenceDerivative = new FiniteDifferenceDerivative(testFunction);
		finiteDifferenceDerivative.setH(0.001);
		NonLinearConjugateGradients nonLinearCG = new NonLinearConjugateGradients();
		nonLinearCG.setEpsilonCG(0.1);
		nonLinearCG.setEpsilonNR(0.1); // 
		nonLinearCG.setImax(1000);
		nonLinearCG.setJmax(5);
		nonLinearCG.nonLinearCG(testFunction, uNew);

		
		System.out.println("u[0] is: "+ uNew[0]);
		System.out.println("u[1] is: "+ uNew[1]);
		
	}
}
