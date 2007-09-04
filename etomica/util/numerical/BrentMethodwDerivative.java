package etomica.util.numerical;

import etomica.util.FunctionDifferentiable;

	/*
	 * 	Section 10.3 Numerical Derivative by Ridder's Method
	 * 	Numerical Recipes in FORTRAN, 2nd Edition, 1992
	 * 
	 *  @ param  double[] x: takes in initial three elements of x; which are the ax, bx, and cx;
	 *  					 in such that ax < bx < cx				
	 *  		  function : from FunctionDifferentiable which is 1-D
	 *  					 f(bx) is less than f(ax) and f(cx)
	 *  			tol    : the fractional precision
	 *  					 the method isolates the minimum to tol using a modification of Brent's method
	 *  					 that uses derivative	         
	 *  
	 *  
	 * 	@ return   xmin   : the minimum double x-value, which is   dbrent[0]
	 * 		      dbrent  : the minimum function value, which is   dbrent[1]
	 * 
	 * 
	 */


public class BrentMethodwDerivative {
	public BrentMethodwDerivative(){
		
	}
	
	public double[] dbrent(double[] x, FunctionDifferentiable function, double tol){
		
		double xmin = 0.0;
		int ITMAX = 100;
		double ZEPS = Math.pow(10, -10);
		
		int iter; 
		boolean ok1, ok2;
		double u;
		double a,b,d1,d2,du,dv,dw,dx; 
		double d = 0.0;
		double e = 0.0;
		double fu,fv,fw,fx,olde,tol1,tol2,u1,u2,v,w,xm;
		
		double ax = x[0];
		double bx = x[1];
		double cx = x[2];
		
		/*
		 * to make sure that a and b are in ascending order
		 */
		a = (ax < cx ? ax: cx);
		b = (ax > cx ? ax: cx);
		
		double axe;
		axe = w = v = bx;
		
		fw= fv = fx = function.f(axe);
		dw= dv = dx = function.df(1, axe);
		
		for (iter=0; iter<ITMAX; iter++){
			xm = 0.5*(a+b);
			tol1 = tol*Math.abs(axe)+ZEPS;
			tol2 = 2.0*tol1;
			
			if(Math.abs(axe-xm) <= (tol2-0.5*(b-a))){
				xmin = axe;
				return new double[]{xmin, fx};
			}
			
			if(Math.abs(e) > tol1){
				d1 = 2.0*(b-a);                     
				d2 = d1;
				
				if(dw != dx) {d1 = (w-axe)*dx/(dx-dw);} 		//Secant Method with one point
				if(dv != dx) {d2 = (v-axe)*dx/(dx-dv);}			// another Second Method
				
				u1 = axe+d1;
				u2 = axe+d2;
				ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
				ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
				
				olde = e;
				e = d;
				
				/*
				 *  Take the acceptable d only; 
				 *  but if double d's are acceptable, only pick the smallest d
				 */
				if(ok1 || ok2){
					if(ok1 &&ok2){
						d = (Math.abs(d1) < Math.abs(d2) ? d1:d2);
					} else if (ok1) {
						d = d1;
					} else {
						d = d2;
					} 
					
					if (Math.abs(d) <= Math.abs(0.5*olde)){
						u = axe + d;
						if (xm-axe >= 0){
							d = Math.abs(tol1);
						} else {
							d = - Math.abs(tol1);
						}
					} else {
						d = 0.5*(e=(dx >= 0.0 ? a-axe : b-axe));
					}
				} else {
					d = 0.5*(e=(dx >= 0.0 ? a-axe : b-axe));
				}
			} else {
				d = 0.5*(e=(dx >= 0.0 ? a-axe : b-axe));
			}
			if (Math.abs(d) >= tol1){
				u = axe+d;
				fu = function.f(u);
			} else {
				if (d >= 0){
					u = axe + Math.abs(tol1);
				} else {
					u = axe - Math.abs(tol1);
				}
				fu = function.f(u);
				
				/*
				 * If the minimum step in the downhill direction takes us uphill,
				 * then we are done.
				 */
				if (fu > fx){
					xmin = axe;
					return new double[] {xmin, fx};
				}
			}
			
			du = function.df(1, u);
			if (fu <= fx){
				if (u >= axe){
					a = axe;
				} else {
					b = axe;
				}
				
				v = w;
				fv = fw;
				dv = dw;
				
				w = axe;
				fw = fx;
				dw = dx;
				
				axe = u;
				fx = fu;
				dx = du;
				
			} else {
				if (u < v){
					a = u;
				} else{
					b = u;
				}
				if (fu <= fw || w == axe){
					
					v = w;
					fv = fw;
					dv = dw;
					
					w = u;
					fw = fu;
					dw = du;
					
				} else if(fu<fv || v==axe || v==w) {
					
					v = u;
					fv = fu;
					dv = du;
					
				}
			}
		}
		
		throw new RuntimeException("Too many iterations in dbrent method!!!");
		//return 0.0;
	}
	
}
