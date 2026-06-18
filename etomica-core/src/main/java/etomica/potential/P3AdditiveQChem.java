/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.util.ParameterBase;

public class P3AdditiveQChem {
	
	public double getU123ADD(P2QChemInterpolated p2, double r12, double x3, double y3) {
		
		double r13 = Math.sqrt(x3*x3 + y3*y3);
		double r23 = Math.sqrt(x3*x3 + (y3-r12)*(y3-r12));
		
		double u12 = p2.u(r12*r12)*kB/JPerHartree;
		double u13 = p2.u(r13*r13)*kB/JPerHartree;
		double u23 = p2.u(r23*r23)*kB/JPerHartree;
    	
		return (u12+u13+u23);
    	
	}

	public static void main(String[] args)  {
		
		P2QChemInterpolated p2 = new P2QChemInterpolated();
		P3AdditiveQChem p3Add = new P3AdditiveQChem();
		
		DampingParams params = new DampingParams();

		if (args.length == 8 ) {
			
			params.r12 = Double.parseDouble(args[0]);
			params.a1 = Integer.parseInt(args[1]);
			params.a2 = Integer.parseInt(args[2]);
			params.RvdwF = Double.parseDouble(args[3]);
			params.basis = Integer.parseInt(args[4]);
			params.fixedRvdw = Boolean.parseBoolean(args[5]);
			params.disp = Boolean.parseBoolean(args[6]);
			params.scf = Boolean.parseBoolean(args[7]);
	    
		} else if (args.length == 1) {
	    	
			params.r12 = Double.parseDouble(args[0]);
			params.a1 = 0;
			params.a2 = 0;
			params.RvdwF = 0;
			params.basis = 0;
			params.fixedRvdw = false;
			params.disp = false;
			params.scf = true;
			
			
		}
		
		int a1 = params.a1;
		int a2 = params.a2;
		double RvdwF = params.RvdwF;
		int basis = params.basis;
		boolean fixedRvdw = params.fixedRvdw;
		boolean disp = params.disp;
		boolean scf = params.scf;
		
		p2.setDampingParams(a1,a2,RvdwF, basis, fixedRvdw);
		p2.setDisp(disp);
		p2.setSCF(scf);
		p2.initialize();
		
		double r12 = params.r12;
   
		
		for (int j=1;j<=100;j++) {

			double x3=j*0.1;


		    for (int k=1;k<=100;k++) {

		    	double y3=k*0.1;
		    	
				double u123Add = p3Add.getU123ADD(p2, r12, x3, y3); //Hartrees
				
				if (u123Add == Double.POSITIVE_INFINITY) {
					System.out.print(String.format("%s\t", "      NaN       ")); 
				} else {
					System.out.print(String.format("%16.15e\t",u123Add));
				}
				 

		    }
		    
		    System.out.println();
		}
		
	   
		
		/*
		double r = 2.3;
		double u = p2.u(r*r);
		System.out.println(r + "  " + u*k/JPerHartree);
    	*/
		//u123Add = p2.u(4*4)+p2.u( 5.65685* 5.65685)+p2.u(4*4);
		
		//System.out.println(u123Add);
		
	}
	
	public static class DampingParams extends ParameterBase {
		
	    protected double r12 = 4;    
	    
	    protected int a1 = 80;	        
	    protected int a2 = 149;   
	    
	    protected int basis = 2;
	    
	    protected double RvdwF = 3.688; 
	    protected boolean fixedRvdw = false; 
	    
	    protected boolean disp = true; 
	    
	    protected boolean scf = false; 
	   
	    
	 

	    
	}
	
	static double kB = 1.3806503e-23; //J/K  
    static double JPerHartree = 4.359744e-18; 

}


