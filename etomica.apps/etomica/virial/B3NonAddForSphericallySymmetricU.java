package etomica.virial;

import etomica.potential.P2QChemInterpolated;
import etomica.potential.P3AdditiveQChem;
import etomica.potential.P3QChem;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.util.ParameterBase;
import etomica.util.SineTransform;


/**
 * 
 * Computes a rough estimate for the classical, non-additive component of B3 for a spherically symmetric potential 
 * 
 * @author kate
 *
 */



public class B3NonAddForSphericallySymmetricU {
	
	public B3NonAddForSphericallySymmetricU() {
		
	}
	
public static void main(String[] args) {

	
	DampingParams params = new DampingParams();
    
   
	if (args.length == 6 ) {
		params.a1 = Integer.parseInt(args[0]);
		params.a2 = Integer.parseInt(args[1]);
		params.Rvdw = Double.parseDouble(args[2]);
		params.basis = Integer.parseInt(args[3]);
		params.fixedRvdw = Boolean.parseBoolean(args[4]);
		params.tempSet = Integer.parseInt(args[5]);
    } 
	
	int a1 = params.a1;
	int a2 = params.a2;
	double Rvdw = params.Rvdw;
	int basis = params.basis;
	boolean fixedRvdw = params.fixedRvdw;
	int tempSet = params.tempSet;
	
	Space space = Space3D.getInstance();
	
	
	
	P2QChemInterpolated p2 = new P2QChemInterpolated(space);
	p2.setDampingParams(a1,a2,Rvdw,basis, fixedRvdw);
	p2.setDisp(true);
	p2.setSCF(true);
	p2.initialize();
	
	P3QChem p3 = new P3QChem();
	
	P3AdditiveQChem p3Add = new P3AdditiveQChem();

	
	double[] temps;
	if (tempSet == 2) {
		temps = new double[] { 100, 133.15, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000 }; // Kelvin
	} else {
		// Primarily temperatures considered by Malejevsky et al:  
		temps = new double[] { 130, 135, 140,145,150,155, 160, 170, 180, 190, 200, 220, 250, 265, 280, 295, 310, 325, 340, 398, 423, 473, 573, 673};
	}
		
		// Some of the temperatures considered by Mas, Lotrich, and Szalewicz (1999):  
		//double[] temps = { 113.15, 133.15, 150.65, 203.15, 323.15, 423.15, 573.16, 673.16, 773.15, 923.15 }; // Kelvin

	
	
	
	System.out.println("T(K)    B3NonAdd (cm6/mol2)");
	
	
	System.out.println();
	
	for (int t=0; t<temps.length; t++) {
		
		double temp = temps[t];
  
		 double B3NonAdd = computeB3NonAdd(p2,p3,p3Add,temp);
			
		System.out.println(temp + "    "  +  B3NonAdd);

	}
	

	}

	public static double computeB3NonAdd(P2QChemInterpolated p2, P3QChem p3, P3AdditiveQChem p3Add, double temp) {
		
		double B3NonAdd = 0;
		for (int i = 25; i<=80; i++) {
			
			double r12 = i*0.1;  //4*Pi*r12*r12 symmetry
			
			for (int j = 1; j<=80; j++) {
				
				double x3 = j*0.1;  // 2*Pi*x symmetry (cylindrical symmetry about y-axis)
				
				for (int k = 0; k<=80; k++) {
					
					double y3 = k*0.1;  
					
					if (y3 >= r12/2) {  // symmetry about r12/2 in xy-plane 
						
						double r13 = Math.sqrt(x3*x3 + y3*y3);
						
						if ((r13 < 8) && (r13 > 2.5) ) {
							
							double r23 = Math.sqrt(x3*x3 - (y3-r12)*(y3-r12));
					
							if ((r23 < 8) && (r23 > 2.5) ) {
						
								//System.out.println(r12+ " " + x3 + "  " + y3);
								
								
								p2.setSCF(false);
								p2.setDisp(true);
								double u123Disp = p3Add.getU123ADD(p2, r12, x3, y3);
								
								
								double u123SCF = p3.getU123SCF(r12, x3, y3); // Hartrees
								
								p2.setSCF(true); // this is set to false to get dispersion energy above
								p2.setDisp(false);
								double u123SCFAdd = p3Add.getU123ADD(p2, r12, x3, y3); // Hartrees
								
								if (u123SCF < -1000) {
									// r12 = 4.7, x3 = 3.1, y3 = 5.7
									// r12 = 5.3, x3 = 3.4, y3 = 4.6
									// r12 = 5.9, x3 = 5.1, y3 = 4.4
									
									u123SCF = Double.NaN;
									
								}
								
								double y3B = y3;
								while (Double.isNaN(u123SCF)) {				
									/*
									y3B = y3B+0.1;
									//System.out.println(y3B + " ");
									if (y3B < y3+0.5) {
										
										u123 = p3.getU123SCF(r12, x3, y3B)+ p3.getU123Disp(p2, r12, x3, y3); // Hartrees
										
										if (u123 < -1000) {
											// r12 = 4.7, x3 = 3.1, y3 = 5.7 : u123SCF = -2635.9 Hartrees --> Infinite e123
											// r12 = 5.3, x3 = 3.4, y3 = 4.6 : u123SCF = -7414.9 Hartrees --> Infinite e123
											// r12 = 5.9, x3 = 5.1, y3 = 4.4 : u123SCF = -3723.6 Hartrees --> Infinite e123
											
											u123 = Double.NaN;
											
										}
										
									} else {
										u123 = u123Add;
									}*/
									u123SCF = u123SCFAdd;
								}
								
								double e123SCF = Math.exp(-u123SCF*JPerH/kB/temp);
								double e123SCFAdd = Math.exp(-u123SCFAdd*JPerH/kB/temp);
								double e123Disp = Math.exp(-u123Disp*JPerH/kB/temp);
								
								double integrand = e123Disp*(e123SCF-e123SCFAdd)*(2*Math.PI*x3*0.1)*(2*0.1)*(4*Math.PI*r12*r12*0.1);
								B3NonAdd = B3NonAdd + integrand;
					
								
								//System.out.println(r12+ " " + x3 + "  " + y3+ "  " + (e123-e123Add));
								
								
							}
						}
						
					}
					
					//
				}
				
				//System.out.println(r12 + " "+ i + " "+ x3 + " "+ j);
			}
			
			//System.out.println(r12);
		}
		
		
		return (-1.0/3.0)*B3NonAdd*0.60221415*0.60221415;
	}
	
	
	
	 public static class DampingParams extends ParameterBase {
		    //TZ
		 	
	    	public int a1 = 80;	        
	        public int a2 = 149;   
	        private double Rvdw = 3.61;
	        private int basis = 2;
	        private boolean fixedRvdw = false;
	        //DZ
	        //public int a1 = 80;	        
	        //public int a2 = 149;   
	       // public int basis = 2;
	        
	        public int tempSet = 3; 
	    }
	 
	 static double kB = 1.3806503e-23; //J/K  
	 static double JPerH = 4.359744e-18; 

}
