/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.util.ParameterBase;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class P3QChem {
	
	public double getU123SCF(double r12, double x3, int k) {
		
		//////////////////////////////////////////////////////
		// Triplet SCF energy (no dispersion)
		//////////////////////////////////////////////////////
		
		double u123SCF = 0;
	
		try{	
			String string;
			
			if ((r12%1.0)==0) {
				string = String.format("%1.0f", r12);
			} else {
				string = String.format("%2.1f", r12);
			}
			
			String file = "/usr/users/kate/ArData/PW86PBE/augccpvdz/matrixSCF"+string+".dat";
			
			//System.out.println(file);
			FileReader fileReader = new FileReader(file);			
			BufferedReader bufReader = new BufferedReader(fileReader);
			String line;
		

		    int j = 0;
		    while ((line = bufReader.readLine()) != null) {
				
		    	j++;
			    	
		    	//System.out.println(j);	
		    	
				double x=j*0.1;
				
				if (x3 == x) {
				
					//int k = (int)Math.round(y3/0.1);
					
					//System.out.println(line);
					String datum = line.split("\t")[k-1];
					//System.out.println(k+"\n"+ datum + y3);
					
					if (datum.length() == 0) {
						u123SCF = Double.NaN;
					} else {
						u123SCF = Double.parseDouble(datum); // Hartrees
					}
					
	
			    }	
		    }
		    
		    bufReader.close();
		    fileReader.close();
		    
		}catch (IOException e){	
			throw new RuntimeException(e);
		}
		
		return u123SCF;
	}

	
	
	public static void main(String[] args)  {
		
		Space space = Space3D.getInstance();
		P2QChemInterpolated p2 = new P2QChemInterpolated();
		
		P3QChem p3 = new P3QChem();
		
		//////////////////////////////////////////////////////
		// Additive dispersion component of triplet energy
		//////////////////////////////////////////////////////
		
		DampingParams params = new DampingParams();
		
		if (args.length == 7 ) {
			
			params.r12 = Double.parseDouble(args[0]);
			params.a1 = Integer.parseInt(args[1]);
			params.a2 = Integer.parseInt(args[2]);
			params.RvdwF = Double.parseDouble(args[3]);
			params.basis = Integer.parseInt(args[4]);
			params.fixedRvdw = Boolean.parseBoolean(args[5]);
			params.disp = Boolean.parseBoolean(args[6]);
	    
		} else if (args.length == 1) {
	    	
			params.r12 = Double.parseDouble(args[0]);
			params.a1 = 0;
			params.a2 = 0;
			params.RvdwF = 0;
			params.basis = 0;
			params.fixedRvdw = false;
			p2.setDampingParams(0,0,0, 0, false);
			params.disp = false;
		}
		
		double r12 = params.r12;
		boolean disp = params.disp;
		if (disp) { 
			int a1 = params.a1;
			int a2 = params.a2;
			double RvdwF = params.RvdwF;
			int basis = params.basis;
			boolean fixedRvdw = params.fixedRvdw;
			
			p2.setDampingParams(a1,a2,RvdwF, basis, fixedRvdw);
			
			p2.setDisp(disp);
			p2.setSCF(false);
			p2.initialize();
		}
		
		P3AdditiveQChem p3Add = new P3AdditiveQChem();
		
		
		
	    for (int j=1;j<=100;j++) {
	    	
	    	double x3 = j*0.1;
	    	
	    	for (int k=1;k<=100;k++) {
		
	    		 double y3 = k*0.1;
	    		 
	    		 double u123 = p3.getU123SCF(r12, x3, k);

	    		 if (disp) {
	    			 u123  = u123 + p3Add.getU123ADD(p2, r12,x3,y3);
	    		 }
	    		 
			    	
		    	if (Double.isNaN(u123)) {
					System.out.print(String.format("%s\t", "      NaN       ")); 
				} else if (Double.isInfinite(u123)) {
					System.out.print(String.format("%s\t", "      NaN       ")); 
				} else{
					System.out.print(String.format("%16.15e\t",u123));
				}
		    	

	    	}
				
			System.out.println();
					
				//System.out.println(j);		
				//if (j > 50) {
					//System.exit(0);
				//}
	    }
		
	}
	
	public static class DampingParams extends ParameterBase {
		
	    protected double r12 = 4;    
	   
	    protected int a1 = 80;	        
	    protected int a2 = 149;   
	    
	    protected int basis = 2;
	    
	    protected double RvdwF = 3.688; 
	    protected boolean fixedRvdw = false; 
	    
	    protected boolean disp = true;
	   
	    
	 

	    
	}
	
	static double kB = 1.3806503e-23; //J/K  
    static double JPerHartree = 4.359744e-18; 

}


