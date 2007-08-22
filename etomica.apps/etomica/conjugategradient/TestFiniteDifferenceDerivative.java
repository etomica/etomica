package etomica.conjugategradient;

import etomica.util.FunctionMultiDimensional;
import etomica.util.numerical.FiniteDifferenceDerivative;


public class TestFiniteDifferenceDerivative implements FunctionMultiDimensional {
	
	public TestFiniteDifferenceDerivative(){
		
	}
	
	public double f(double[] u){
        return 2*u[0] + u[0]*u[1]*u[1]*u[1];
		//return Math.cos(u[0]);
	}
	
	public int getDimension() {
	    return 1;   
    }


	public static void main(String args[]){
		
		
		double[] u = new double[]{1,20};//rad
		TestFiniteDifferenceDerivative test = new TestFiniteDifferenceDerivative();
		
		FiniteDifferenceDerivative differentiate = new FiniteDifferenceDerivative(test);
		differentiate.setH(0.001);
		differentiate.setHOptimizer(true);
		System.out.println("cosine of 45 is "+test.f(u));
		System.out.println("Differentiated Result is: "+ differentiate.df(new int[] {1,2}, u));
        System.out.println("analytic: " + 3*u[1]*u[1]);
		// -sin(45) should be -0.850903524 rad
		
		
	}
}