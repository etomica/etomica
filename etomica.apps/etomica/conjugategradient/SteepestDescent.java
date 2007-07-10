package etomica.conjugategradient;

import Jama.Matrix;

public class SteepestDescent {
	protected int imax;
	protected Matrix A; 
	protected Matrix B, X, R, Q; 
	protected double delta, delta0;
	
	public SteepestDescent(int iter, Matrix a, Matrix x, Matrix b){
		this.imax = iter;
		this.B = b;
		this.A = a;
		this.X = x;
	}
	
	public Matrix SteepestDescentAlgorithm(){
		//double [] val = new double []{};
		int i = 0;
		
		/*
		 * an error tolerance, 0 <= epsilon < 1
		 * as default, let epsilon equals to 0.5
		 */
		double epsilon = 0.5;
		
		// r = b - A*x
		R = B;
		R.minusEquals(A.times(X));
		
		// delta = r^T*r
		delta = ((R.transpose()).times(R)).trace();
		delta0 = delta;
		
		
		while (i<imax && delta > epsilon*epsilon*delta ){
			Q = A.times(R);
			double alpha = delta/(((R.transpose()).times(Q)).trace());
			
			X.plusEquals(R.timesEquals(alpha));
			
			/*
			 * every square root of number of iterations, exact residual is computed
			 * As default, set every 50 iterations
			*/
			if (i%50 == 0){
				R = B;
				R.minusEquals(A.times(X));
			} else {
				R.minusEquals(Q.timesEquals(alpha));
			}
			
			delta = ((R.transpose()).times(R)).trace();
			
			i++;
		}
		
		System.out.println(X.det());
		return X;
	}

	public static void main (String [] args){
		
	      double[][] valsA = {{1.,3.},{0.375,1.}};
	      Matrix aValue = new Matrix(valsA);
	      
	      double[][] valsB = {{2.},{1}};
	      Matrix bValue = new Matrix(valsB);
	      
	      double[][] valsX = {{6.},{-1.}};
	      Matrix xValue = new Matrix(valsX);
	      
	      SteepestDescent steepestDescent = new SteepestDescent(10000, aValue, xValue, bValue);
	      steepestDescent.SteepestDescentAlgorithm();
	    //  System.out.println("xValue: " + xValue);
		
	}
}

