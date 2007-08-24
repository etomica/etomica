package etomica.conjugategradient;

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
	
	
	public NonLinearConjugateGradients(FiniteDifferenceDerivative finiteDifferenceDerivative){
		this.finiteDifferenceDerivative = finiteDifferenceDerivative;
		
	}
	
	public void NonLinearCG(int imax, int jmax, double epsilon, double[] u){
		
		
		/*
		 *  imax is a maximum number of CG iterations
		 *  jmax is a maximum number of Newton-Raphson iterations
		 *  
		 */
		
		int i=0;
		int k=0;
		int n;
		
		double deltaNew = 0;
		int coordinateDim = u.length;
		int[] dAssign = new int[coordinateDim];
		double[] uDerivative = new double[coordinateDim];
		
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
			uDerivative[diff1Counter] = finiteDifferenceDerivative.df(dAssign, u);
		}
		
		
		double[] r = new double[coordinateDim]; 
		double[] d = new double[coordinateDim];
		
		for (n=0; n<coordinateDim; n++){
			r[n] = - uDerivative[n];
			d[n] = r[n];
			
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
			int[][] d2Assign = new int[coordinateDim][coordinateDim];
			double[][] u2Derivative = new double[coordinateDim][coordinateDim];
			//////////////////////
			
			for (int diff2Counter =0; diff2Counter< coordinateDim; diff2Counter++){
				for (int diff1Counter =0; diff1Counter< coordinateDim; diff1Counter++){
					/*
					 *  To assign d to differentiate over all the dimensions
					 */
					
					
					
					for (int dim =0; dim< coordinateDim; dim++){
						if (diff2Counter==dim){
							d2Assign[dim][dim] = 1;
						} else{
							d2Assign[dim][dim] = 0;
						}
					}
					
					u2Derivative[diff2Counter][diff1Counter] = finiteDifferenceDerivative.df(dAssign, uDerivative);
				}
				
			}
			
			
			///////////////////////
			
			double[] d_DoublePrime = new double[coordinateDim];
			
			for(n=0; n<coordinateDim; n++){
				
				for(int m=0; m<coordinateDim; m++){
					d_DoublePrime[n] += d[m]*u2Derivative[m][n];
				}
			}
			
			for(n=0; n<coordinateDim; n++){
				
				deltad += d[n]*d[n];
				
				alpha_num += - uDerivative[n]*d[n];
				alpha_denom += d_DoublePrime[n]*d[n];
			}
			
			double alpha = alpha_num /alpha_denom;
			
			for(n=0; n<coordinateDim; n++){
				u[n] = u[n] + alpha*d[n];
			}
			j++;
			
			double alpha2_deltad = alpha*alpha*deltad;
			double epsilon2 = epsilon*epsilon;
			
			while(j<jmax && alpha2_deltad > epsilon2){
				double deltaOld=0;
				
				for(n=0; n<coordinateDim; n++){
					r[n] = - uDerivative[n];
					deltaOld = deltaNew;
					
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
	
}
