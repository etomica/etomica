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
	 */
	
	public void mnbrak(double[] x, FunctionDifferentiable function){
		double GOLD = 1.618034; 
		double GLIMIT = 100.0;
		double TINY = Double.MIN_VALUE;
		double denom = 0.0;
		
		double ax = x[0];
		double bx = x[1];
		double cx = x[2];
		
		double dum=0, u;
		double fu, q, r, ulimit;
		
		double fa = function.f(ax);
		double fb = function.f(bx);
		
		/*
		 *  Switching the roles of a and b so that able to go
		 *  downhill in the direction from a to b 
		 */
		if (fb >= fa){
			
			dum = ax;		
			ax = bx;
			bx = dum;
			
			dum = fb;
			fb = fa;
			fa = dum;
		}
		
		/*
		 * Initial guess for cx
		 * GOLD is 1 plus the fraction of the so-called Golden Mean or Golden Section 
		 * 	(0.38197 and 0.61893 respectively, if considered from both different ends)
		 */
		cx = bx + GOLD*(bx-ax);
		double fc = function.f(cx);
		
		while (fb > fc){
			r = (bx - ax)*(fb - fc);
			q = (bx - cx)*(fb - fa);
			
			if ((q-r) >=0){
				denom = 2*Math.abs(Math.max(Math.abs(q-r), TINY));
			} else {
				denom = -2*Math.abs(Math.max(Math.abs(q-r), TINY));
			}
			
			/*
			 * Computing u by parabolic extrapolation from ax, bx, and cx
			 * TINY is used to avoid division by zero
			 */
			u = bx - ((bx-cx)*q - (bx-ax)*r)/denom;
			ulimit = bx + GLIMIT*(cx-bx);
			
			if((bx-u)*(u-cx) > 0.0){ // If Parabolic u is between b and c
				fu = function.f(u);
				if (fu < fc){        // Got a minimum between b and c
					ax = bx;
					bx = u;
					fa = fb;
					fb = fu;
					break;
				} else if (fu > fb){ // Got a minimum between a and u
					cx = u;
					fc = fu;
					break;
				}
				
				/*
				 * When Parabolic fit was not used
				 * Use default magnification
				 */
				
				u = cx + GOLD*(cx-bx);
				fu = function.f(u);
			
			} else if ((cx-u)*(u-ulimit) > 0.0){        // If Parabolic fit is between c and its allowed limit
				fu = function.f(u);
				if (fu < fc){
					
					bx = cx;
					cx = u;
					u  = cx+GOLD*(cx-bx);
					
					fb = fc;
					fc = fu;
					fu = function.f(u);
				}
				
			} else if ((u-ulimit)*(ulimit-cx) >= 0.0){   // Limitting parabolic u to the maximum allowed value
				u = ulimit;
				fu = function.f(u);
				
			} else {									 // Rejecting parabolic u, use default magnification									
				u = cx+ GOLD*(cx-bx);
				fu = function.f(u);
			}
			
			/*
			 *  Eliminating oldest point and continue
			 */
			
			ax = bx;
			bx = cx;
			cx = u;
			
			fa = fb;
			fb = fc;
			fc = fu;
			
		}
		
		x[0] = ax;
		x[1] = bx;
		x[2] = cx;
	}
	
}
