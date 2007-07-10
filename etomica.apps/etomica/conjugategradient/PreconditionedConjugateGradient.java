package etomica.conjugategradient;

import Jama.Matrix;

public class PreconditionedConjugateGradient {

	protected int imax;
	protected Matrix A, M, D, S; 
	protected Matrix B, X, R, Q; 

	
	public PreconditionedConjugateGradient(Matrix A, Matrix x, Matrix b, Matrix m, int iter){
		
		this.A = A;
		this.X = x;
		this.B = b;
		this.M = m;
		this.imax = iter;
		
	}
	
	
   /**  Preconditioned Conjugate Gradient solution for
    * 	Ax = b
   @param A     Matrix (square, symmetric, and positive-definite/indefinite)
   @param x		Vector (initial solution guess)
   @param b		Vector
   @param M     Symmetric, positive-definite matrix that approximates A
   
   @return x 	Solution (five significant figures)
   
   */
	public Matrix PreconditionedConjugateGradientAlgorithm(){
			
		/*
		 * an error tolerance, 0 <= epsilon < 1
		 * as default, let epsilon equal to 0.000001
		 */
		double epsilon = 0.0000001;
		
		double deltaNew, deltaOld;
		double delta0;
		double alpha;
		double e2d0;
		
		int i = 0;
		
		// r = b - A*x
		R = B.copy();
		R.minusEquals(A.times(X));
			
		// d = M^(-1) * r
		D = M.inverse().times(R);
		
		deltaNew = ((R.transpose()).times(D)).trace();
		delta0 = deltaNew;
		
		e2d0 = epsilon*epsilon*delta0;
		
		while (i < imax && deltaNew > e2d0 ){
		
			Q = A.times(D);
			
			alpha =  deltaNew / ((D.transpose()).times(Q)).trace();
					
			X.plusEquals(D.times(alpha));
			
			//Check for exact residual computation or fast recursive formula.  As default, exact set every 50 iterations.
			if (i%50 == 0){
				R = B.copy();
				R.minusEquals(A.times(X));
			} 
			
			else if (deltaNew <= e2d0){
				R = B.copy();
				R.minusEquals(A.times(X));
			}
			else {
				R.minusEquals(Q.times(alpha));
			}
			
			S = M.inverse().times(R);
			deltaOld = deltaNew;
			deltaNew = (R.transpose().times(S)).trace();
			
			double beta = deltaNew /deltaOld;
			
			D = S.copy();
			D.plusEquals(D.times(beta));
			
			i++;
			
		}
		return X;
	}

	public static void main (String [] args){
		
	      double[][] valsA = {{1.,3.},{.375,1.}};
	      Matrix A = new Matrix(valsA);
	      
	      double[][] valsB = {{2.},{1.}};
	      Matrix b = new Matrix(valsB);
	      
	      double[][] valsX = {{-7.},{1.}};
	      Matrix x = new Matrix(valsX);
	      
	      double[][] valsM = {{2.,1.},{1.,2.}};
	      Matrix m = new Matrix(valsM);
	      
	      //A.print(10,3);
	      //b.print(10,3);
	      //x.print(10,3);
	      
	      PreconditionedConjugateGradient preconditionedConjugateGradient = new PreconditionedConjugateGradient(A, x, b, m, 1000);
	      
	      x = preconditionedConjugateGradient.PreconditionedConjugateGradientAlgorithm();
	      
	      x.print(10, 4);
		
	}
}
