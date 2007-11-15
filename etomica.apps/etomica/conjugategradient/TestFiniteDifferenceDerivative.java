package etomica.conjugategradient;

import etomica.util.FunctionMultiDimensional;
import etomica.util.numerical.FiniteDifferenceDerivative;


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
		differentiate.setHOptimizer(true);
		differentiate.setNtab(10);
		
		System.out.println("Differentiated Result is: "+ differentiate.df(new int[] {2, 0}, u));
		System.out.println("Value for function is: "+test.f(u));
		
		
	}
}