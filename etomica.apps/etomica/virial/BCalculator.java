package etomica.virial;

import etomica.math.SpecialFunctions;
import etomica.virial.cluster.ClusterBz;

/**
 * @author kofke
 *
 * Calculates the density-expansion virial coefficients given values for the
 * activity-expansion coefficients.
 */

//need to check that coefficient of b2^8 is correct for B9

public class BCalculator {

	/**
	 * Constructor for BCalculator.
	 */
	public BCalculator() {
		super();
	}
	
	public double[] getB() {
		return B;
	}
	
	public void setBz(double[] bz) {
		int n = bz.length;
		B = new double[n];
		
		for(int m=2; m<n+2; m++) {
			ClusterBz.IntegerSet iSet = new ClusterBz.IntegerSet(m);
			B[m-2] = 0.0;
			do {
				if(!acceptableSet(iSet.set)) continue;
				double coeff = coefficient(iSet.set);
				double prod = 1.0;
				System.out.println("m, c, set: "+m+" "+coeff+" "+iSet.toString());
				for(int j=2; j<=m; j++) prod *= power(bz[j-2],iSet.set[j-1]);
				B[m-2] += coeff*prod;
			} while(iSet.advanceFull());
			System.out.println();
		}
	}
	
	private static boolean acceptableSet(int[] set) {
		int sum1 = 0;
		int sum2 = 0;
		final int m = set.length;
		for(int j=1; j<=m; j++) {
			sum1 += set[j-1];
			sum2 += j*set[j-1];
		} 
		return (sum1 == m-1) && (sum2 == 2*(m-1));
	}
	
	private static int coefficient(int[] set) {
		final int m = set.length;
		int sign = (((m-set[0]-1) % 2)==0) ? +1 : -1; //(-1)^(m-k1-1)
		int coeff = sign * (m-1) * SpecialFunctions.factorial(2*m-set[0]-3);
		for(int j=2; j<=m; j++) {//product of kj terms
			int kj = set[j-1];
			if(kj == 0) continue;
			int kjFact = SpecialFunctions.factorial(kj);
			int jkj = power(j,kj);
			coeff *= jkj;
			coeff /= kjFact;
		}
		coeff /= SpecialFunctions.factorial(m);
		return coeff;	
	}
	
	//returns k^n
	private static int power(int k, int n) {
		if(n < 0) throw new IllegalArgumentException();
		int prod = 1;
		for(int i=n; i>0; i--) prod *= k;
		return prod;
	}
	
	private static double power(double a, int n) {
		if(n < 0) throw new IllegalArgumentException();
		double prod = 1.0;
		for(int i=n; i>0; i--) prod *= a;
		return prod;		
	}
	private double[] B = new double[0];

	public static void main(String[] args) {
		BCalculator calc = new BCalculator();
		double b2 = 5.31575;
		double[] bRatios = new double[] { 1.0    ,2.03392  ,1.55684   ,1.0171};
//		double b2 = 1.3145;
//		double[] bRatios = new double[] { 1.0    ,2.76782  ,4.35115   ,165.072};
//		double b2 = 3.31774;
//		double[] bRatios = new double[] { 1.0    ,2.1245  ,1.7378   ,1.23096};
//		double[] bError = new double[]  { 0.0    ,0.002  ,0.002  ,0.01};
		
		double b2HS = -2.*Math.PI/3.;
		double b3RatioHS = -8.0/(5./8. - 4.);
		double[] bRatiosHS = new double[] {1.0, 2.3703703703703703704, 
												2.2509954663107567552, //b4
												1.8892976717277971446, //b5
												1.466374, //b6
												1.0823, //b7
												0.757};//b8
		double[] bError = new double[] {1.0, b3RatioHS, 
												0.003, //b4
												0.004, //b5
												0.007, //b6
												0.01, //b7
												0.007};//b8
		bRatios = bRatiosHS; b2 = b2HS;
		double[] bz = new double[bRatios.length]; //bk = bz[k-2]
		bz[0] = b2;
//		for(int i=1; i<bz.length; i++) bz[i] = (1.0/bRatios[i])*power(2.0*b2,i+1);
		for(int i=1; i<bz.length; i++) bz[i] = (1.0/bRatios[i])*power(2,i+1); bz[0] = 1.0;
		calc.setBz(bz);
		double[] B = calc.getB();
		for(int i=0; i<B.length; i++) {
			System.out.println("B"+(i+2)+": "+B[i]+"  "+B[i]/power(B[0],i+1));
		}
	}
}
