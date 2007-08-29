package etomica.conjugategradient;

import etomica.util.FunctionMultiDimensional;
import etomica.util.numerical.FiniteDifferenceDerivative;

public class TestNonLinearConjugateGradients implements FunctionMultiDimensional{
	
	public TestNonLinearConjugateGradients(){
		
	}
	
	public double f(double[] u){
		return u[0]*u[0] + u[1]*u[1];
	}
	
	public int getDimension(){
		return 1;
	}
	
	public static void main(String args[]){
		double[] uNew = new double[] {.5, .5};
		TestNonLinearConjugateGradients testFunction = new TestNonLinearConjugateGradients();
		
		FiniteDifferenceDerivative finiteDifferenceDerivative = new FiniteDifferenceDerivative(testFunction);
		finiteDifferenceDerivative.setH(0.001);
		finiteDifferenceDerivative.setHOptimizer(true);
		NonLinearConjugateGradients nonLinearCG = new NonLinearConjugateGradients(finiteDifferenceDerivative);
		nonLinearCG.setEpsilonCG(0.1);
		nonLinearCG.setEpsilonNR(0.71); // 
		nonLinearCG.setImax(1000);
		nonLinearCG.setJmax(10);
		nonLinearCG.NonLinearCG(testFunction, uNew);

		
		System.out.println("u[0] is: "+ uNew[0]);
		System.out.println("u[1] is: "+ uNew[1]);
		
	}
}
