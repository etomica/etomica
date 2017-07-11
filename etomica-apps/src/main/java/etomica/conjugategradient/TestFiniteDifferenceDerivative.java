/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.conjugategradient;

import etomica.math.function.FunctionMultiDimensional;
import etomica.math.numerical.FiniteDifferenceDerivative;


public class TestFiniteDifferenceDerivative implements FunctionMultiDimensional {
	
	public TestFiniteDifferenceDerivative(){
		
	}
	
	public double f(double[] u){
        return u[0]*u[0] + u[1]*u[1];
		//return Math.sin(u[0]);
	}
	
	public int getDimension() {
	    return 1;   
    }


	public static void main(String args[]){
		
		
		double[] u = new double[]{1, 1};
		TestFiniteDifferenceDerivative test = new TestFiniteDifferenceDerivative();
		
		FiniteDifferenceDerivative differentiate = new FiniteDifferenceDerivative(test);
		differentiate.setH(0.0001);
		
		System.out.println("Differentiated Result is: "+ differentiate.df(new int[] {2, 0}, u));
		System.out.println("Value for function is: "+test.f(u));
		
		
	}
}
