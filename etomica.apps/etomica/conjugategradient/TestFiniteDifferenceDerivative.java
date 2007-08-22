package etomica.conjugategradient;

import etomica.util.FunctionMultiDimensional;


public class TestFiniteDifferenceDerivative implements FunctionMultiDimensional {
	
	public TestFiniteDifferenceDerivative(){
		
	}
	
	public double function(double[] u){
		
		return Math.cos(u[0]);
	}
	
	


	public static void main(String args[]){
		
		
		double[] u = new double[]{45};//rad
		TestFiniteDifferenceDerivative test = new TestFiniteDifferenceDerivative();
		
		FiniteDifferenceDerivativeWOBox differentiate = new FiniteDifferenceDerivativeWOBox(test);
		differentiate.setH(0.001);
		differentiate.setHOptimizer(false);
		System.out.println("cosine of 45 is "+test.function(u));
		System.out.println("Differentiated Result is: "+ differentiate.dfdx(u)[0]);
		// -sin(45) should be -0.850903524 rad
		
		
	}
}