package etomica.conjugategradient;

import java.io.OutputStream;

import Jama.Matrix;

public class TestMatrix {
	public static void main(String[] args){
	      double[][] vals = {{1.,2.,3},{4.,5.,6.},{7.,8.,10.}};
	      Matrix A = new Matrix(vals);
	      Matrix b = Matrix.random(3,1);
	      Matrix x = A.solve(b);
	      Matrix r = A.times(x).minus(b);
	      double rnorm = r.normInf();
	      
	      
	}
}
