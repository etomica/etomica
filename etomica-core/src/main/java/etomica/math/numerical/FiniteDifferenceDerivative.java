/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.numerical;

import etomica.math.function.FunctionMultiDimensional;
import etomica.math.function.FunctionMultiDimensionalDifferentiable;

public class FiniteDifferenceDerivative implements FunctionMultiDimensionalDifferentiable {
	
	protected FunctionMultiDimensional fFunction;
	protected double h;
	
	public FiniteDifferenceDerivative(FunctionMultiDimensional fFunction){

		this.fFunction = fFunction;
		h = 0.00001;
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
        int[] dCopy = d.clone();
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
        
        double uSave = u[index];

        u[index] = uSave + h;
        double fPlus = df(dCopy, u);
        u[index] = uSave - h;
        double fMinus= df(dCopy, u);
            
        double a = (fPlus - fMinus)/(2.0*h);

        //System.out.println(" a[0][0] is: "+a[0][0]);
        u[index] = uSave;
        return a;
	}

	public double getH() {
		return h;
	}

	public void setH(double h) {
		this.h = h;
	}

}
