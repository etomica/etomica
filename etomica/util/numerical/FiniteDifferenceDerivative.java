package etomica.util.numerical;

import etomica.util.FunctionMultiDimensional;
import etomica.util.FunctionMultiDimensionalDifferentiable;

public class FiniteDifferenceDerivative implements FunctionMultiDimensionalDifferentiable {
	
	/*
	 * Section 5.7 Numerical Derivative by Ridder's Method
	 * Numerical Recipes in FORTRAN, 2nd Edition, 1992
	 * 
	 * @author Tai Tan
	 */
		
	protected FunctionMultiDimensional fFunction;
	protected double h;
	protected boolean hOptimizer;
	protected int ntab;
	
	public FiniteDifferenceDerivative(FunctionMultiDimensional fFunction){

		this.fFunction = fFunction;
		h = 0.00001;
		hOptimizer = false;
		ntab = 10;
	}
	
	public double f(double[] u){
		return fFunction.f(u);
	}
    
    public int getDimension() {
        return fFunction.getDimension();
    }

	public double df(int[] d, double[] u) {
		
		/*
		 * 
		 */
		
        if(u.length != d.length) {
            throw new IllegalArgumentException("d and u must be the same length");
        }
        
        int index = -1;
        int[] dCopy = (int[])d.clone();
        for(int i=0; i<d.length; i++) {
            if(d[i] != 0) {
                index = i;
                dCopy[i] = d[i] - 1;
                break;
            }
        }
        if(index == -1) {
            return fFunction.f(u);
        }
        
		double con = 1.4;
		double con2 = con*con;
		double big = Double.MAX_VALUE;
		double safe = 2.0;
		
		double errt, fac, hh;
		double[][] a = new double[ntab][ntab];
		
		hh = h;

        double uSave = u[index];

        u[index] = uSave + hh;
        double fPlus = df(dCopy, u);
        u[index] = uSave - hh;
        double fMinus= df(dCopy, u);
            
        a[0][0] = (fPlus - fMinus)/(2.0*hh);

        //System.out.println(" a[0][0] is: "+a[0][0]);
		if (!hOptimizer) {
            u[index] = uSave;
			return a[0][0];
		}
			
		double err = big;
        double dfdx = Double.NaN;
		
		for(int i=1; i<ntab; i++){
			hh = hh /con;
            u[index] = uSave + hh;
            fPlus = df(dCopy, u);
            u[index] = uSave - hh;
            fMinus= df(dCopy, u);
			
            a[0][i] = (fPlus - fMinus)/(2.0*hh);
			//System.out.println(" a[0]["+i+"] is: "+a[0][i]);
			fac = con2;
			
			for(int j=1; j<i; j++){
				a[j][i] = (a[j-1][i]*fac - a[j-1][i-1])/(fac-1);
				//System.out.println(" a["+j+"]["+i+"] is: "+a[j][i]);
				fac = con2*fac;
				errt = Math.max(Math.abs(a[j][i]-a[j-1][i]), Math.abs(a[j][i]-a[j-1][i-1]));
				//System.out.println("errt is: "+errt);
				
				if (errt <= err){
					err = errt;
					dfdx = a[j][i];
					//System.out.println("in errt<= err, dfdx is: "+a[j][i]);
				}
			}
			
			if (Math.abs(a[i][i]-a[i-1][i-1]) >= safe*err){
				break;
			}
		}
		u[index] = uSave;
		return dfdx;
	}

	public double getH() {
		return h;
	}

	public void setH(double h) {
		this.h = h;
	}

	public boolean isHOptimizer() {
		return hOptimizer;
	}

	public void setHOptimizer(boolean optimizer) {
		hOptimizer = optimizer;
	}

	public int getNtab() {
		return ntab;
	}

	public void setNtab(int ntab) {
		this.ntab = ntab;
	}
}
