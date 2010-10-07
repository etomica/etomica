package etomica.util.numerical;


/**
 * 
 * This class is coded for specific nonlinear model that decribes the Equation:
 *  f(x) = polynomialA + exp(-Kx)*polynomialB
 * 
 * M and N are the input parameters that determine the order of polynomialA
 *  and polynomialB respectively
 *  
 * polynomialA is given as: 
 *                M-1
 *  polynomialA = sum [ a_m * x^m ]
 *                m=0
 * 
 * for example, when M=3; polynomialA = a_0 + a_1 * x + a_2 * x^2
 * 
 * 
 * Note: polynomialB has same form as polynomialA 
 *               
 * @author Tai Boon Tan
 *
 */
public class FittingFunctionNonLinear{
	public FittingFunctionNonLinear(int M, int N){
		this.M = M;
		this.N = N;
	}

	public double f(double[] a, double x) {
		double sumPolyA = 0.0;
		double sumPolyB = 0.0;
		double expKx = Math.exp(-a[M]*x);

		//polynomialA
		double xMultiply = 1;
		for(int i=0; i<(M-1); i++){
			sumPolyA += a[i]*xMultiply;
			xMultiply *= x;
		}
		
		//polynomialB
		xMultiply = 1;
		for(int i=(M+1); i<(M+N+1); i++){
			sumPolyB += a[i]*xMultiply;
			xMultiply *= x;
		}
		
		return sumPolyA + expKx*sumPolyB;
		
	}

	/*
	 * d is derivatives
	 * a is array of parameters
	 * x is the x-value
	 */
	public double df(int d, double[] a, double x) {
		double dsumPolyA = 0.0;
		double dsumPolyB = 0.0;
		double dexpKx = Math.exp(-a[M]*x);
		double xMultiply;
		
		if(d==M){
			dexpKx *= -x;
		}
		
		//polynomialA
		xMultiply = 1;
		for(int i=0; i<(M-1); i++){
			if(i==d){
				dsumPolyA += xMultiply;
			} else {
				dsumPolyA += a[i]*xMultiply;
			}
			xMultiply *= x;
		}
		
		//polynomialB
		xMultiply = 1;
		for(int i=(M+1); i<(M+N+1); i++){
			if(i==d){
				dsumPolyB += xMultiply;
			} else {
				dsumPolyB += a[i]*xMultiply;
			}
			xMultiply *= x;
		}
		
		return dsumPolyA + dexpKx*dsumPolyB;

	}
	
	protected int M, N;
}
