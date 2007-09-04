package etomica.util.numerical;

import etomica.util.FunctionMultiDimensionalDifferentiable;

public class TestConjugateGradientsMultiDimensional implements FunctionMultiDimensionalDifferentiable{
	
	public TestConjugateGradientsMultiDimensional(){
		
	}
	
	public double f(double[] u){
		//return u[0]*u[0] + u[1]*u[1];                              // x ^2 + y ^2
		//return (u[0]+1)*(u[0]+1) + (u[1]+1)*(u[1]+1);            // (x + 1) ^2 + (y + 1) ^2 
		return u[0]*u[0]+(u[0]+u[1])*(u[0]+u[1]);                // x ^2 + (x + y) ^2
		//return Math.exp(u[0])+Math.exp(-u[0]) + u[1]*u[1];         // e^(x) + e^(-x) + y^2
	}
	
	public double df(int[] d, double[] u){
		//		 x ^2 + y ^2
//		if (d[0]==0){
//			if (d[1]==0){
//				return f(u);
//			} else if (d[1]==1){
//				return 2*u[1];
//			} else if (d[1]==2){
//				return 2;
//			} 
//		} else if (d[0]==1){
//			return(d[1]==0)? 2*u[0]: 0.0;
//		} else if (d[0]==2){
//			return(d[1]==0)? 2: 0.0;
//		}
		
		
		//		 (x + 1) ^2 + (y + 1) ^2  
//		if (d[0]==0){
//			if (d[1]==0){
//				return f(u);
//			} else if (d[1]==1){
//				return 2*(u[1]+1);
//			} else if (d[1]==2){
//				return 2;
//			} 
//	} else if (d[0]==1){
//		return(d[1]==0)? 2*(u[0]+1): 0.0;
//	} else if (d[0]==2){
//		return(d[1]==0)? 2: 0.0;
//	}
		
		
		//		 x ^2 + (x + y) ^2
		if (d[0]==0){
			if (d[1]==0){
				return f(u);
			} else if (d[1]==1){
				return 2*u[0]+2*u[1];
			} else if (d[1]==2){
				return 2;
			} 
	} else if (d[0]==1){
		return(d[1]==0)? 4*u[0]+2*u[1]: 2.0;
	} else if (d[0]==2){
		return(d[1]==0)? 4: 0.0;
	}
		
		//		 e^(x) + e^(-x) + y^2
//		if (d[0]==0){
//			if (d[1]==0){
//				return f(u);
//			} else if (d[1]==1){
//				return 2*u[1];       // dy
//			} else if (d[1]==2){
//				return 2;            // dy dy
//			} 
//		} else if (d[0]==1){
//			return(d[1]==0)? Math.exp(u[0]) - Math.exp(-u[0]): 0.0;  // dx : dx dy
//		} else if (d[0]==2){
//			return(d[1]==0)? Math.exp(u[0]) + Math.exp(-u[0]): 0.0;       // dx dx
//		}
		
		
		return 0;
		
	}
	
	public int getDimension(){
		return this.getDimension();
	}
	
	public static void main(String args[]){
		double[] p = new double[] {1, 1};
		TestConjugateGradientsMultiDimensional testFunction = new TestConjugateGradientsMultiDimensional();
		
		ConjugateGradientMultiDimensional conjugateGradient = new ConjugateGradientMultiDimensional();
		conjugateGradient.conjugateGradient(p, 0.00001, testFunction);
		
		System.out.println("Minimum Function value is: "+ conjugateGradient.getFunctionMinimimumValue());
		
		for(int i=0; i<p.length; i++){
			System.out.println("Minimum value of p["+i+"] is: "+ conjugateGradient.getMinimumCoordinates()[i]);
		}
		
		System.out.println("Number of iteration is: "+ conjugateGradient.getNumIterations());
		
		
		//System.out.println("u[0] is: "+ uNew[0]);
		//System.out.println("u[1] is: "+ uNew[1]);
		
	}
}
