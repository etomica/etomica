package etomica.util.numerical;

/*
 * unmodified Besselfuntion or modified first type
 *  
 * author Weisong Lin
 */

public class BesselFunction {

	public static double J(int alpha, double x) {
		return doCalc(false, alpha, x);
	}

	public static double I(int alpha, double x) {
		return doCalc(true, alpha, x);
	}
	
	public static double doCalc(boolean modified ,int alpha,double x){//true means modified fisrt type and false means unmodified
		double sum = 0;
		double am = 0;
		double factorial1 = 1;
		double factorial2 = 1;

		if(modified){
			double s1 = x*x/4;
			double s2 = Math.pow(x/2, alpha);
			for(int m=0;m<1E10;m++){

				factorial1 *= (m==0)? 1:m;
				factorial2 *= ((m+alpha)==0)?1:(m+alpha);
				am = 	s2/(factorial1*factorial2) ;
				
//				System.out.println("f1 = " + factorial1+ " f2=" + factorial2 + " a"+ m+" = " + am );
//				System.out.println("a"+ m+" = " + am );
				
				if((am/sum)<1E-15){
//					System.out.println("ratio=" + (am/sum));
//					System.out.println("sum= " + sum );
					break;
				}
				sum +=  am;
				s2 *= s1;
			}
			
		}else{
			int s = 1;
			for(int m=0;m<1E10;m++){
				factorial1 *= (m==0)? 1:m;
				factorial2 *= ((m+alpha)==0)?1:(m+alpha);
				am = 	s*Math.pow(x/2, 2*m+alpha)/(factorial1*factorial2) ;//TODO simpliy this one!
				if((am/sum)<1E-15){
					break;
				}
				sum +=  am;
				s = -s;
			}

		}
		
		return sum;
	}
	
	
}
