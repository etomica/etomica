package etomica.conjugategradient;

import Jama.Matrix;

/**
 * @author taitan
 * @author msellers
 *
 */

public class SteepestDescent {

	protected int imax;
	protected Matrix A; 
	protected Matrix B, X, R, Q; 

	
	public SteepestDescent(Matrix A, Matrix x, Matrix b, int iter){
		this.imax = iter;
		this.B = b;
		this.A = A;
		this.X = x;
	}
	
	
   /**  Steepest Descent solution for
    * 	Ax = b
   @param A     Matrix (square, symmetric, and positive-definite/indefinite)
   @param x		Vector (initial solution guess)
   @param b		Vector
   @return x 	Solution (five significant figures)
   */
	public Matrix SteepestDescentAlgorithm(){
			
		/*
		 * an error tolerance, 0 <= epsilon < 1
		 * as default, let epsilon equal to 0.0001
		 */
		double epsilon = 0.000001;
		
		double delta;
		double delta0;
		double alpha;
		double e2d0;
		
		int i = 0;
		
		// r = b - A*x
		R = B.copy();
		R.minusEquals(A.times(X));
			
		// delta = r^T*r
		delta = ((R.transpose()).times(R)).trace();
		
		delta0 = delta;
		e2d0 = epsilon*epsilon*delta0;
		
		while (i < imax && delta > e2d0 ){
		
			Q = A.times(R);
			
			alpha =  delta / ((R.transpose()).times(Q)).trace();
					
			X.plusEquals(R.times(alpha));
			
			System.out.println("delta = "+delta);
			System.out.println("e2d0 = "+e2d0);
			
			//Check for exact residual computation or fast recursive formula.  As default, exact set every 50 iterations.
			if (i%50 == 0){
				R = B.copy();
				R.minusEquals(A.times(X));
			} 
			
			else if (delta <= e2d0){
				R = B.copy();
				R.minusEquals(A.times(X));
			}
			else {
				R.minusEquals(Q.times(alpha));
			}
			
			delta = ((R.transpose()).times(R)).trace();
			
			i++;
			
		}
		
		return X;
	}

	public static void main (String [] args){
		
	      double[][] valsA = {{3.,2.},{2.,6.}};
	      Matrix A = new Matrix(valsA);
	      
	      double[][] valsB = {{2.},{-8.}};
	      Matrix b = new Matrix(valsB);
	      
	      double[][] valsX = {{-8.},{8.}};
	      Matrix x = new Matrix(valsX);
	      
	      A.print(10,3);
	      b.print(10,3);
	      x.print(10,3);
	      
	      SteepestDescent steepestDescent = new SteepestDescent(A, x, b, 150);
	      
	      x = steepestDescent.SteepestDescentAlgorithm();
	      
	      x.print(10, 4);
		
	}
}

