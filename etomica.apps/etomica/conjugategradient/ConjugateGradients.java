package etomica.conjugategradient;

import Jama.Matrix;

/**
 * @author taitan
 * @author msellers
 *
 */

public class ConjugateGradients {

	protected int imax;
	protected Matrix A; 
	protected Matrix B, X, R, Q; 
	
	
	
	public ConjugateGradients(Matrix A, Matrix x, Matrix b, int iter){
		this.imax = iter;
		this.B = b;
		this.A = A;
		this.X = x;
	}
	
	
   /**  Conjugate Gradients solution for
    * 	Ax = b
   @param A     Matrix (square, symmetric, and positive-definite/indefinite)
   @param x		Vector (initial guess for unknown solution)
   @param b		Vector
   @return x 	Vector (solution)
   */
	public Matrix ConjugateGradientsAlgorithm(){
		
		/*
		 * an error tolerance, 0 <= epsilon < 1
		 * as default, let epsilon equal to 0.0001
		 */
		double epsilon = 0.000001;
		
		double deltaNEW;
		double deltaOLD;
		double delta0;
		double alpha;
		double beta;
		double e2d0;
		Matrix D;
		
		int i = 0;
		
		// r = b - A*x
		R = B.copy();
		R.minusEquals(A.times(X));
			
		D = R.copy();
		
		// deltaNEW = r^T*r
		deltaNEW = ((R.transpose()).times(R)).trace();
		
		delta0 = deltaNEW;
		e2d0 = epsilon*epsilon*delta0;
		
		while (i < imax && deltaNEW > e2d0 ){
		
			Q = A.times(D);
			
			alpha =  deltaNEW / ((D.transpose()).times(Q)).trace();
					
			X.plusEquals(D.times(alpha));
			
			//Check for exact residual computation or fast recursive formula.  As default, exact set every 50 iterations.
			if (i%50 == 0){
				R = B.copy();
				R.minusEquals(A.times(X));
			} 
			
			else if (deltaNEW <= e2d0){
				R = B.copy();
				R.minusEquals(A.times(X));
			}
			else {
				R.minusEquals(Q.times(alpha));
			}
			
			deltaOLD = deltaNEW;
			deltaNEW = ((R.transpose()).times(R)).trace();
			
			beta = deltaNEW/deltaOLD;
			
			D = (R.copy()).plusEquals(D.times(beta));
			
			i++;
			
		}
		
		return X;
	}

	public static void main (String [] args){
		
	      double[][] valsA = {{3.,2.},{2.,6.}};
	      Matrix A = new Matrix(valsA);
	      
	      double[][] valsB = {{2.},{-8.}};
	      Matrix b = new Matrix(valsB);
	      
	      double[][] valsX = {{-26.},{253.}};
	      Matrix x = new Matrix(valsX);
	      
	      A.print(10,3);
	      b.print(10,3);
	      x.print(10,3);
	      
	      ConjugateGradients conjugateGradient = new ConjugateGradients(A, x, b, 150);
	      
	      x = conjugateGradient.ConjugateGradientsAlgorithm();
	      
	      x.print(10, 4);
	}
}
