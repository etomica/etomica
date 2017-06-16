/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.math.numerical.AkimaSpline;

/**
 * 
 * Soft-Sphere Solid Equation of State
 *  for n = 12, 9 and 6
 * 
 * Reference: Tan, Schultz and Kofke, J Chem. Phys.,133(2010)
 * 			   Appendix: Free Energy Formulas
 * 
 * @author Tai Boon Tan
 *
 */
public class SoftSphereSolidEOS {
	
	public SoftSphereSolidEOS(){
		akimaSpline = new AkimaSpline();
	}
		
	public double getUlat0(int id, int n){
		return ulat0[id][n/3-2];
	}
	
	public double getVib0(int id, int n){
		return vib0[id][n/3-2];
	}
	
	public double[] getCoeffAc(int id , int n){
		return c[id][n/3-2];
	}
	
	public double[] getAcQuantity(int id, int n, double temp, double rho){
		double bAc = 0.0;
		double dAcdrho = 0.0;
		double dbAcdbeta = 0.0;
		
		double[] c = getCoeffAc(id, n);
		double upsilon = getUpsilon(n, temp, rho);
		double ups = 1; 
		
		for(int i=1; i<=7; i++){
			ups *= upsilon;
			bAc += c[i]*ups;
			
			dAcdrho += temp*c[i]*ups*(-n*i/(3.0*rho));
			dbAcdbeta += c[i]*ups*(-i*temp);
		}
		
		double[] values = new double[]{bAc, dAcdrho, dbAcdbeta};
		
		return values;
	}
	
	public double[] getAllQuantity(int id, int n, double temp, double rho){
		double upsilon = getUpsilon(n, temp, rho);
		double ulat0 = getUlat0(id, n);
		double Vib0 = getVib0(id, n);
		double[] AcValues = getAcQuantity(id, n, temp, rho); 
		
		double bUlat = ulat0/upsilon; 			            // Equation (A2)
		double bAVib = Vib0 - 3.0/2.0*Math.log(upsilon) 
		                + Math.log(rho*sigma*sigma*sigma);  // Equation (A3)
		double bAc   = AcValues[0];			 				// Equation ( 9)
		double bA    = bUlat + bAVib + bAc; 				// Equation (A1) 
		
		double p = rho*rho*((n/3.0)*(temp/(upsilon*rho))*ulat0  	
				   + AcValues[1]
				   + (n*temp)/(2*rho) + 1/rho);
		
		double u = temp/upsilon*ulat0
				   + AcValues[1]
				   + (3.0/2.0)*temp;
		
		// retur betaA, Z (P/rho*kT), betaU
		double[] quantity = new double[]{bA, p/(rho*temp), u/temp};
		return quantity;
	}
	
	public double getUpsilon(int n, double temp, double rho){
		return (temp/epsilon)*Math.pow((rho*sigma*sigma*sigma), -n/3.0);
	}
	
	public double getPressure(int id, int n, double temp, double rho){
		double upsilon = getUpsilon(n, temp, rho);
		double ulat0 = getUlat0(id, n);
		double[] AcValues = getAcQuantity(id, n, temp, rho); 
		
		double p = rho*rho*((n/3.0)*(temp/(upsilon*rho))*ulat0  	
				   + AcValues[1]
				   + (n*temp)/(2*rho) + 1/rho);
		
		return p;
	}
	
	public double findDensityatP(double desiredP, int id, int n, double temp){
		
		double[] p = new double[20];
		double[] density = new double[20];
		double rhoMin, rhoMax;
		
		if(n==12){rhoMin=1.1; rhoMax=3.0;}
		else if(n==9){rhoMin=1.2; rhoMax=3.0;}
		else {rhoMin=2.0; rhoMax=3.0;};
		
		int interval = density.length;
		
		for(int i=0; i<interval; i++){
			density[i] = rhoMin+i*(rhoMax-rhoMin)/interval; 
//			System.out.println(i+" "+density[i]);
		}
		
        for(int i=0; i<density.length; i++){
        	p[i] = getPressure(id, n, temp, density[i]);
//        	System.out.println(i+" "+p[i]);
        }
        
        akimaSpline.setInputData(p, density);
        double[] rho =  akimaSpline.doInterpolation(new double[]{desiredP});
        double del = 0.0003;

        while (del > 1e-10) {
          p = new double[5];
          double[] newRho = new double[]{rho[0]-2*del, rho[0]-del, rho[0], rho[0]+del, rho[0]+2*del};
       
          for(int i=0; i<p.length; i++){
        	  p[i] = getPressure(id, n, temp, density[i]);
          	  
          }
          akimaSpline.setInputData(p, newRho);
          rho =  akimaSpline.doInterpolation(new double[]{desiredP});
          del *= 0.1;
        }
		return rho[0];
	}
	
	public double findTransP(int idA, int idB, int n, double temp, int initP){
		
		double P0 = initP;
		double rhoA = findDensityatP(P0, idA, n, temp);
		double rhoB = findDensityatP(P0, idB, n, temp);
		double[] qA = getAllQuantity(idA, n, temp, rhoA);
		double[] qB = getAllQuantity(idB, n, temp, rhoB);
		
		double bGA0 = qA[0] + qA[1];
		double bGB0 = qB[0] + qB[1];
		double deltabG0 = bGB0 - bGA0;
	
		double P1 = initP + 40;
		rhoA = findDensityatP(P1, idA, n, temp);
		rhoB = findDensityatP(P1, idB, n, temp);
		qA = getAllQuantity(idA, n, temp, rhoA);
		qB = getAllQuantity(idB, n, temp, rhoB);
		
		double bGA1 = qA[0] + qA[1];
		double bGB1 = qB[0] + qB[1];
		double deltabG1 = bGB1 - bGA1;
	
//		System.out.println("deltabG0: " + deltabG0);
//		System.out.println("deltabG1: " + deltabG1);
//		System.exit(1);
		double P=Double.NaN;
		double rhoANew = Double.NaN;
		double rhoBNew = Double.NaN;
		double deltabGNew = 1000;
		int i =0;
		while (Math.abs(deltabGNew) > tolerance){
			
			// Interpolation to determine P at deltaG = 0.0
			P = P0 + ((P1-P0)/(deltabG1-deltabG0))*(0.0 - deltabG0);
			System.out.println("P: " + P);
			
			rhoANew = findDensityatP(P, idA, n, temp);
			qA = getAllQuantity(idA, n, temp, rhoANew);
			double bGANew = qA[0] + qA[1];
			
			
			rhoBNew = findDensityatP(P, idB, n, temp);
			qB = getAllQuantity(idB, n, temp, rhoBNew);
			double bGBNew = qB[0] + qB[1];
	
			deltabGNew = bGBNew - bGANew;
			
//			System.out.println("deltabGNew: " + deltabGNew);
//			System.exit(1);
			if(deltabGNew <= 0.0){
				deltabG0 = deltabGNew;
				P0 = P;
			} else {
				deltabG1 = deltabGNew;
				P1 = P;
			}	
			
			if(i==2000){
				throw new RuntimeException("CANNOT find transition properties!!!!!!");
			}
			
			++i;
		}
		
		return P;
	}
	
	public static void main (String[] args){
		double temp = 1.0;
		
		SoftSphereSolidEOS ssA = new SoftSphereSolidEOS();
		System.out.println("--------- SOFT SPHERE SOLID EQUATION OF STATE FOR N=12, 9 AND 6 -----------\n");
		// Find beta.A/N, compressibility factor, Z and beta.U/N at desired density
		if(true){
			int n = 6;
			double rho = 2.32; 
			
			double[] qFCC = ssA.getAllQuantity(FCC, n, temp, rho);
			double[] qBCC = ssA.getAllQuantity(BCC, n, temp, rho);
			System.out.println("**** betaA, Z and betaU for n="+n+ " at rho="+rho + " ****");
			System.out.println("[FCC] betaA: " + qFCC[0] + "; Z: " + qFCC[1] + "; betaU: "+ qFCC[2]);
			System.out.println("[BCC] betaA: " + qBCC[0] + "; Z: " + qBCC[1] + "; betaU: "+ qBCC[2]+"\n");
//			return;
		}
		
		// Find Gibbs Free Energy at constant pressure
		if(true){
			int n = 6;
			double desiredP = 114;
			
			System.out.println("****  betaG for n="+n+ " at P="+desiredP + "  ****");
			double rhoFCC = ssA.findDensityatP(desiredP, FCC, n, temp);
			double rhoBCC = ssA.findDensityatP(desiredP, BCC, n, temp);
			double[] qFCC = ssA.getAllQuantity(FCC, n, temp, rhoFCC);
			double[] qBCC = ssA.getAllQuantity(BCC, n, temp, rhoBCC);

			System.out.println("[FCC] betaG: "+ (qFCC[0]+qFCC[1]) + " at density: "+rhoFCC);
			System.out.println("[BCC] betaG: "+ (qBCC[0]+qBCC[1]) + " at density: "+rhoBCC+"\n");
//			return;
		}
		
		// Find transition pressure for fcc and bcc (n=6)
		if(false){
			double transP = ssA.findTransP(FCC, BCC, 6, 1.0, 90);
			System.out.println("transition pressure: " + transP);
		}
		
	}
	
	protected AkimaSpline akimaSpline;
	protected double tolerance = 1e-10;
	
    public final static int FCC = 0;
    public final static int HCP = 1;
    public final static int BCC = 2;
    protected double epsilon = 1.0;
    protected double sigma = 1.0;
    
    protected final double[][] ulat0 = new double[][]{{3.613475382,2.208391122,1.516485025},{3.613718378,2.208528128,1.516536721},{3.630711328}};
    protected final double[][] vib0 = new double[][]{{2.69819834,3.48217524,3.88845245},{2.70269996,3.48590058,3.89158992},{2.56473110}};
    protected final double[][][] c = new double[][][]{
    		{{0.0, -0.024270, -0.385254, -1.989540, 24.366184, -208.856577, 860.698810,-1469.431853},
    		 {0.0,  0.249526, -0.245724, -0.500979,  4.258323,  -16.542027,  30.811960,  -23.208566},
    		 {0.0,  0.464148, -0.423821,  0.466322, -0.526094,    0.148125,   0.401893,   -0.384095}},
    		
    		{{0.0, -0.023687, -0.351963, -2.934478, 37.070278, -280.028273, 953.287096, -1140.820226},
    		 {0.0, 0.250099,  -0.282838,  0.007132,  0.922672,   -5.218832,  10.890988,    -7.369508},
    		 {0.0, 0.462065,  -0.413562,  0.429175, -0.461716,    0.111491,   0.368595,    -0.352244}},
    		
    		{{0.0, 0.283913, -2.421458, 14.558777, -77.105689, 229.711416, -391.161702, 335.089915}},
    			 
    };
    
}
