package etomica.util.numerical;


/**
 * 
 * 
 * 
 * @author taitan
 *
 */
public class FittingFunctionNonLinear{
	public FittingFunctionNonLinear(){
		
	}

	public double f(double[] a, double x) {
		return a[0]*Math.exp(-a[1]*x);
		
	}

	/*
	 * d is derivatives
	 * a is array of parameters
	 * x is the x-value
	 */
	public double df(int d, double[] a, double x) {
		if(d==0){
			return Math.exp(-a[1]*x);
			
		} else if (d==1){
			return -a[0]*x*Math.exp(-a[1]*x);
			
		} else {
			return Double.NaN;
			
		}
	}
	
}
