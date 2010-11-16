package etomica.util.numerical;

import etomica.util.FunctionDifferentiable;

public class BisectionMethodMinimumBracket {

	public BisectionMethodMinimumBracket(){
		
	}
	
	/*
	 * 	Section 10.2 Numerical Derivative by Ridder's Method
	 * 	Numerical Recipes in FORTRAN, 2nd Edition, 1992
	 * 
	 * Given a function and distinct initial points ax and bx
	 * 	The method searches in the downhill direction
	 * 	returns new points ax, bx and cx that bracket a minimum of the function
	 * 
	 * @param double[] x: takes in initial three elements of x; which are the ax, bx, and cx;
	 * 					  cx in this case should be zero.								
	 * 		  function  : from FunctionDifferentiable which is 1-D
	 * 
	 * Author: Tai Boon Tan
	 */
	
	public double[] mnbrak(double[] x, FunctionDifferentiable function){
		double GOLD = 1.618034; 
		double GLIMIT = 100.0;
		double TINY = Double.MIN_VALUE;
		double denom = 0.0;
		
		double[] newX = x;
		
		newX[0] = x[0];
		newX[1] = x[1];
		newX[2] = x[2]; 
		
		double dum=0, u;
		double fu, q, r, ulimit;
		
		double fa = function.f(newX[0]);
		double fb = function.f(newX[1]);
		

		/*
		 *  Switching the roles of a and b so that able to go
		 *  downhill in the direction from a to b 
		 */
		if (fb > fa){
			
			dum = newX[0];		
			newX[0] = newX[1];
			newX[1] = dum;
			
			dum = fb;
			fb = fa;
			fa = dum;
		}
		
		/*
		 * Initial guess for newX[2]
		 * GOLD is 1 plus the fraction of the so-called Golden Mean or Golden Section 
		 * 	(0.38197 and 0.61893 respectively, if considered from both different ends)
		 */
		newX[2] = newX[1] + GOLD*(newX[1]-newX[0]);
		double fc = function.f(newX[2]);

		while (fb >= fc){
			r = (newX[1] - newX[0])*(fb - fc);
			q = (newX[1] - newX[2])*(fb - fa);
			
			if ((q-r) >=0){
				denom = 2*Math.abs(Math.max(Math.abs(q-r), TINY));
			} else {
				denom = -2*Math.abs(Math.max(Math.abs(q-r), TINY));
			}
			
			/*
			 * Computing u by parabolic extrapolation from newX[0], newX[1], and newX[2]
			 * TINY is used to avoid division by zero
			 */
			u = newX[1] - ((newX[1]-newX[2])*q - (newX[1]-newX[0])*r)/denom;
			ulimit = newX[1] + GLIMIT*(newX[2]-newX[1]);
			
			if((newX[1]-u)*(u-newX[2]) > 0.0){ // If Parabolic u is between b and c
				fu = function.f(u);
				if (fu < fc){        // Got a minimum between b and c
					newX[0] = newX[1];
					newX[1] = u;
					fa = fb;
					fb = fu;
					return newX;
				} else if (fu > fb){ // Got a minimum between a and u
					newX[2] = u;
					fc = fu;
					return newX;
				}
				
				/*
				 * When Parabolic fit was not used
				 * Use default magnification
				 */
				
				u = newX[2] + GOLD*(newX[2]-newX[1]);
				fu = function.f(u);
			
			} else if ((newX[2]-u)*(u-ulimit) > 0.0){        // If Parabolic fit is between c and its allowed limit
				fu = function.f(u);
				if (fu < fc){
					
					newX[1] = newX[2];
					newX[2] = u;
					u  = newX[2]+GOLD*(newX[2]-newX[1]);
					
					fb = fc;
					fc = fu;
					fu = function.f(u);
				}
				
			} else if ((u-ulimit)*(ulimit-newX[2]) >= 0.0){   // Limitting parabolic u to the maximum allowed value
				u = ulimit;
				fu = function.f(u);
				
			} else {									 // Rejecting parabolic u, use default magnification									
				u = newX[2]+ GOLD*(newX[2]-newX[1]);
				fu = function.f(u);
			}
			
			/*
			 *  Eliminating oldest point and continue
			 */
			
			newX[0] = newX[1];
			newX[1] = newX[2];
			newX[2] = u;
			
			fa = fb;
			fb = fc;
			fc = fu;
			
		}
		
		x[0] = newX[0];
		x[1] = newX[1];
		x[2] = newX[2];

		return newX;
	}
	
}
