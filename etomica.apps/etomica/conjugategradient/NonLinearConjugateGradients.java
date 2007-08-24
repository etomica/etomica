package etomica.conjugategradient;

import etomica.util.FunctionMultiDimensionalDifferentiable;
import etomica.util.numerical.FiniteDifferenceDerivative;


public class NonLinearConjugateGradients {

	/*
	 *   Nonlinear Conjugate Gradients with Newton-Raphson and Fletcher-Reeves
	 *  "An Introduction to the Conjugate Gradient Method w/o the Agonizing Pain"
	 *                      Jonathan Richard Shewchuk
	 *                           August 4, 1994
	 *                           
	 *  @author Tai Tan
	 */
	
	protected FiniteDifferenceDerivative finiteDifferenceDerivative;
	protected int imax;
	protected int jmax;
	protected double epsilon;
	
	public NonLinearConjugateGradients(FiniteDifferenceDerivative finiteDifferenceDerivative){
		this.finiteDifferenceDerivative = finiteDifferenceDerivative;
		imax = 100;
		jmax = 100;
		epsilon = 0.00001;
	}
	
	public void NonLinearCG(FunctionMultiDimensionalDifferentiable function, double[] u){
		
		/*
		 *  imax is a maximum number of CG iterations
		 *  jmax is a maximum number of Newton-Raphson iterations
		 *  
		 */
		
		int i=0;
		int k=0;
		int n;
		int coordinateDim = u.length;
		double deltaNew = 0;
		
		int[] dAssign = new int[coordinateDim];
		double[] uDerivative = new double[coordinateDim];
		
		finiteDifferenceDerivative = new FiniteDifferenceDerivative(function);
		
		/*
		 * To obtain first derivative of a function
		 */
		for (int diff1Counter=0; diff1Counter< coordinateDim; diff1Counter++){
			/*
			 *  To assign d to differentiate over all the dimensions
			 */
			for (int dim=0; dim< coordinateDim; dim++){
				if (diff1Counter==dim){
					dAssign[dim] = 1;
				} else{
					dAssign[dim] = 0;
				}
			}
			uDerivative[diff1Counter] = - finiteDifferenceDerivative.df(dAssign, u);
		}
		
		
		double[] r = uDerivative.clone(); 
		double[] d = r.clone();
		
		for (n=0; n<coordinateDim; n++){	
			System.out.println("r["+n + "] is: " + r[n]);
			deltaNew += r[n]*r[n];
		}
		
		double delta0 = deltaNew;
		double epsilon2_delta0 = epsilon*epsilon*delta0;
		
		System.out.println("NonlinearCG before WHILE loop...");
		System.out.println("DeltaNew: "+ deltaNew);
		System.out.println("epsilon: "+ epsilon);
		System.out.println("epsilon2_delta0: "+ epsilon2_delta0);
		
		while(i<imax && deltaNew > epsilon2_delta0){
			int j=0;
			double deltad = 0;
			double alpha_num = 0;
			double alpha_denom = 0;
			
			/*
			 * To obtain second derivative of a function
			 */
			double[][] u2Derivative = new double[coordinateDim][coordinateDim];
		
			for (int diff2Counter =0; diff2Counter< coordinateDim; diff2Counter++){
				for (int diff1Counter =0; diff1Counter< coordinateDim; diff1Counter++){
					
					/*
					 *  To assign d to differentiate twice over all the dimensions
					 */
					for (int dim =0; dim< coordinateDim; dim++){
						if (diff1Counter==dim){
							dAssign[dim] = 1;
						} else{
							dAssign[dim] = 0;
						}
						
						++ dAssign[diff2Counter];  // determine the next element of the derivative
					}
					
					u2Derivative[diff2Counter][diff1Counter] = finiteDifferenceDerivative.df(dAssign, u);
				}
			}
			
			double[] d_DoublePrime = new double[coordinateDim];
			
			/*
			 * Computing Matrix-Multiplication of 
			 * Matrix (1 X CoordinateDim) and Matrix (CoordinateDim X CoordinateDim)
			 */
			for(n=0; n<coordinateDim; n++){
				
				for(int m=0; m<coordinateDim; m++){
					d_DoublePrime[n] += d[m]*u2Derivative[m][n];
				}
			}
			
			for(n=0; n<coordinateDim; n++){
				
				deltad += d[n]*d[n];
				
				alpha_num += uDerivative[n]*d[n];
				alpha_denom += d_DoublePrime[n]*d[n];
			}
			
			double alpha = - alpha_num /alpha_denom;
			
			for(n=0; n<coordinateDim; n++){
				u[n] = u[n] + alpha*d[n];
			}
			j++;
			
			double alpha2_deltad = alpha*alpha*deltad;
			double epsilon2 = epsilon*epsilon;
			
			while(j<jmax && alpha2_deltad > epsilon2){
				
				r = uDerivative.clone();
				double deltaOld = deltaNew;
				
				for(n=0; n<coordinateDim; n++){
					deltaNew += r[n]*r[n];
				}
				
				double beta = deltaNew/ deltaOld;
				
				for(n=0; n<coordinateDim; n++){
					d[n] = r[n] + beta*d[n];
				}
				
				k++;
			}
			
			double rTd = 0;
			
			for(n=0; n<coordinateDim; n++){
				rTd += r[n]*d[n];
			}

			if(k==n || rTd <= 0){
				for(n=0; n<coordinateDim; n++){
					d[n] = r[n];
				}
				k=0;
			}
			
			System.out.println("NonlinearCG within WHILE loop...");
			i++;
		}
	}

	public int getImax() {
		return imax;
	}

	public void setImax(int imax) {
		this.imax = imax;
	}

	public int getJmax() {
		return jmax;
	}

	public void setJmax(int jmax) {
		this.jmax = jmax;
	}

	public double getEpsilon() {
		return epsilon;
	}

	public void setEpsilon(double epsilon) {
		this.epsilon = epsilon;
	}
	
}
