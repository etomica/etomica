/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.normalmode.ArrayReader1D;
import etomica.units.Pascal;
import etomica.util.Constants;
import etomica.math.numerical.AkimaSpline;
import etomica.math.numerical.PolynomialFit;

/**
 * Class that determines the transition properties for alpha and beta nitrogen
 * 
 * @author taitan
 *
 */


public class FindTransitionTemperature {
	
	public FindTransitionTemperature(double[] density, double desiredPressure, double tolerance, String fname){
		akimaSpline = new AkimaSpline();
		this.density = density;
		this.desiredPressure = desiredPressure;
		this.tolerance = tolerance;
		
		alphaCoeffs = new double[density.length][6];
		betaCoeffs  = new double[density.length][5];
		
		for(int i=0; i<density.length; i++){
			alphaCoeffs[i] = ArrayReader1D.getFromFile(fname+"alphaCoeffd"+density[i])[0];
			betaCoeffs[i] = ArrayReader1D.getFromFile(fname+"betaCoeffd"+density[i])[0];
		}
		
	}
	
	public double[] AandPAlpha(double rho, double temperature){
		int numDensityToUse = checkNumDensityToUse(temperature); 
		double[] a = new double[numDensityToUse];
		double[] densityLimit = new double[numDensityToUse];
		
		for(int i=0; i<a.length; i++){
			a[i] = AAlphaFunction(i, temperature);
			densityLimit[i] = density[i];
//			System.out.println("a[" + i+"]: " + a[i]);
//			System.out.println(density[i]+" " + a[i]);

		}
//		System.exit(1);
		/*
		 * Smoothen out the curve by fitting a[i] to 
		 *  3rd order polynomial 
		 */
		
		double[] r = PolynomialFit.doFit(3, densityLimit, a);
		
		double A = r[0] + r[1]*rho + r[2]*rho*rho + r[3]*rho*rho*rho;
		double P = r[1] + 2*r[2]*rho + 3*r[3]*rho*rho;
		double[] AandP = new double[]{A, P*rho*rho}; 

//		System.out.println("rho: " + rho);
//		System.out.println("temperature: " + temperature);
//		System.out.println(rho + " "+aAlpha[0]+" "+pAlpha[0]);
//		System.exit(1);
		return  AandP;
	}
	
	public double[] AandPBeta(double rho, double temperature){
		double[] a = new double[density.length];
//		System.out.println("\n"+rho +" ***A**");
		for(int i=0; i<a.length; i++){
			a[i] = ABetaFunction(i, temperature);
//			System.out.println("a[" + i+"]: " + a[i]);
//			System.out.println(density[i]+" " + a[i]);
		}
//		System.exit(1);
		/*
		 * Smoothen out the curve by fitting a[i] to 
		 *  4th order polynomial 
		 */
		double[] r = PolynomialFit.doFit(4, density, a);

		double A = r[0] + r[1]*rho + r[2]*rho*rho + r[3]*rho*rho*rho + r[4]*rho*rho*rho*rho;
		double P = r[1] + 2*r[2]*rho + 3*r[3]*rho*rho + 4*r[4]*rho*rho*rho;
		double[] AandP = new double[]{A, P*rho*rho}; 
		
		return  AandP;
	}
	
	public double UAlpha(double TransRho, double TransT){
		double[] betaA = new double[41];
		double[] invTemp = new double[41];
		
		for(int i=0; i<betaA.length; i++){
			double temp = 49.0-i*((49.0-36.0)/(betaA.length-1));
			invTemp[i] = 1.0/temp;
			betaA[i] = AandPAlpha(TransRho, temp)[0]*invTemp[i];
//			System.out.println(invTemp[i] + " " + betaA[i]);
		}
		
		double invT = 1.0/TransT;
		akimaSpline.setInputData(invTemp, betaA);
		double[] u = akimaSpline.doInterpolationDy(new double[]{invT});
		
//		for(int i=0; i<u.length; i++){
//			System.out.println(invTemp[i]+" " + u[i]);
//		}
		
//		System.out.println("alpha energy at " + invT  + " " + u[0]);
		return  u[0];
	}
	
	
	public double UBeta(double TransRho, double TransT){
		double[] betaA = new double[41];
		double[] invTemp = new double[41];
		
		for(int i=0; i<betaA.length; i++){
			double temp = 49.0-i*((49.0-36.0)/(betaA.length-1));
			invTemp[i] = 1.0/temp;
			betaA[i] = AandPBeta(TransRho, temp)[0]*invTemp[i];
//			System.out.println(invTemp[i] + " " + betaA[i]);
		}
		
		double invT = 1.0/TransT;
		akimaSpline.setInputData(invTemp, betaA);
		double[] u = akimaSpline.doInterpolationDy(new double[]{invT});
			
//		System.out.println("beta energy at " + invT  + " " + u[0]);
		return  u[0];
	}
	
	public double AAlphaFunction(int i, double temperature){
		double T=temperature;
		double A = 0.0;
		
		double[] coeff = alphaCoeffs[i];
		A = coeff[0];
		for (int k=1; k<coeff.length; k++){
			A += coeff[k]*T;
			T *=temperature;
		}

		return A;
		
	}
	
	public double ABetaFunction(int i, double temperature){
		double T=temperature;
		double A = 0.0;
		
		double[] coeff = betaCoeffs[i];
		A = coeff[0];
		for (int k=1; k<coeff.length; k++){
			A += coeff[k]*T;
			T *=temperature;
		}

		return A;
				
	}
		
	public double findDensityAlpha(double temperature){
	    int numDensityToUse = checkNumDensityToUse(temperature); 
		double [] p = new double[numDensityToUse];
	    
        for(int i=0; i<p.length; i++){
        	p[i] = AandPAlpha(density[i], temperature)[1];
//        	System.out.println(density[i]+" "  + p[i]);
//        	System.out.println("p["+i+"]: " + p[i]);
        }
//        System.exit(1);
        akimaSpline.setInputData(p, density);
        double[] rho =  akimaSpline.doInterpolation(new double[]{desiredPressure});
        double del = 0.0003;
        
        while (del > 1e-10) {
          p = new double[5];
          double[] newRho = new double[]{rho[0]-2*del, rho[0]-del, rho[0], rho[0]+del, rho[0]+2*del};

          for(int i=0; i<p.length; i++){
        	  p[i] = AandPAlpha(newRho[i], temperature)[1];
//        	  System.out.println("p["+i+"]: " + p[i]);
          }

          akimaSpline.setInputData(p, newRho);
          rho =  akimaSpline.doInterpolation(new double[]{desiredPressure});
//          System.out.println("rho: " + rho[0]);
//          System.exit(1);
          del *= 0.1;
        }
//        System.out.println("density: " + rho[0]);
		return rho[0];
	}
	
	public double findDensityBeta(double temperature){
		
	    double [] p = new double[density.length];

        for(int i=0; i<p.length; i++){
        	p[i] = AandPBeta(density[i], temperature)[1];
//        	System.out.println("p["+i+"]: " + p[i]);
//        	 System.out.println(density[i] +" " + p[i]);
        }
//        System.exit(1);
        akimaSpline.setInputData(p, density);
        double[] rho =  akimaSpline.doInterpolation(new double[]{desiredPressure});
        double del = 0.0003;

        while (del > 1e-10) {
          p = new double[5];
          double[] newRho = new double[]{rho[0]-2*del, rho[0]-del, rho[0], rho[0]+del, rho[0]+2*del};
       
          for(int i=0; i<p.length; i++){
        	  p[i] = AandPBeta(newRho[i], temperature)[1];
//        	  System.out.println("p["+i+"]: " + p[i]);
//          	  System.out.println(density[i] +" " + p[i]);
          	  
          }
//          System.exit(1);
          akimaSpline.setInputData(p, newRho);
          rho =  akimaSpline.doInterpolation(new double[]{desiredPressure});
          del *= 0.1;
        }
//        System.out.println(temperature + " density: " + rho[0]);
		return rho[0];
	}
	
	public double findTransTemperature(double T0){
		
		//Alpha-phase Gibbs free energy
		// working with per particle quantity) 
		// GAlpha = G/N or newAlphaA = A/N
		double T1;
		
		
//		for (int temp=49; temp<50; temp++){	
//			double AlphaRho0 = findDensityAlpha(temp);
//			double[] AandPAlpha0 = AandPAlpha(AlphaRho0, temp);
//			double AlphaG0 = AandPAlpha0[0] + AandPAlpha0[1]/AlphaRho0; 
//			
////			double BetaRho0 = findDensityBeta(temp);
////			double[] AandPBeta0 = AandPBeta(BetaRho0, temp);
////			double BetaG0 = AandPBeta0[0] + AandPBeta0[1]/BetaRho0; 
//		
////			double deltaG0 = BetaG0 - AlphaG0;
//			System.out.println(temp + " "+AandPAlpha0[0] + " " + AandPAlpha0[1] + " "+ AlphaRho0+ " "+AlphaG0);
////			System.out.println(AandPBeta0[0] + " " + AandPBeta0[1] + " "+ BetaRho0+ " "+BetaG0);
////			System.out.println(temp + " deltaG0: " + deltaG0);
//		}
//		System.exit(1);
		
		
		double AlphaRho0 = findDensityAlpha(T0);
//		System.out.println("AlphaRho0: " + AlphaRho0);
//		System.exit(1);
		double[] AandPAlpha0 = AandPAlpha(AlphaRho0, T0);
		double AlphaG0 = AandPAlpha0[0] + AandPAlpha0[1]/AlphaRho0; 
		
		double BetaRho0 = findDensityBeta(T0);
		double[] AandPBeta0 = AandPBeta(BetaRho0, T0);
		double BetaG0 = AandPBeta0[0] + AandPBeta0[1]/BetaRho0; 
	
		double deltaG0 = BetaG0 - AlphaG0;
//		System.out.println(AandPAlpha0[0] + " " + AandPAlpha0[1] + " "+ AlphaRho0+ " "+AlphaG0);
//		System.out.println(AandPBeta0[0] + " " + AandPBeta0[1] + " "+ BetaRho0+ " "+BetaG0);
//		System.out.println("deltaG0: " + deltaG0);
//		System.exit(1);
		
		if(deltaG0 >=0.0){
			T1 = 44;
			
		} else {
			T1 = 36;
		}
//		System.out.println("T1: " + T1);
		double AlphaRho1 = findDensityAlpha(T1);
		double[] AandPAlpha1 = AandPAlpha(AlphaRho1, T1);
		double AlphaG1 = AandPAlpha1[0] + AandPAlpha1[1]/AlphaRho1; 
		
		double BetaRho1 = findDensityBeta(T1);
		double[] AandPBeta1 = AandPBeta(BetaRho1, T1);
		double BetaG1 = AandPBeta1[0] + AandPBeta1[1]/BetaRho1; 
		
		
//		System.exit(1);
		double deltaG1 = BetaG1 - AlphaG1;
//		System.out.println(AandPAlpha1[0] + " " + AandPAlpha1[1] + " "+ AlphaRho1+ " "+AlphaG1);
//		System.out.println(AandPBeta1[0] + " " + AandPBeta1[1] + " "+ BetaRho1+ " "+BetaG1);
//		System.out.println("deltaG1: " + deltaG1);
		
		double T=Double.NaN;
		double AlphaRhoNew = Double.NaN;
		double BetaRhoNew = Double.NaN;
		double deltaGNew = 1000;
		int i =0;
		while (Math.abs(deltaGNew) > tolerance){
			
			// Interpolation to determine T at deltaG = 0.0
			T = T0 + ((T1-T0)/(deltaG1-deltaG0))*(0.0 - deltaG0);
//			System.out.println("T: " + T);
			AlphaRhoNew = findDensityAlpha(T);
			double[] AandPAlphaNew = AandPAlpha(AlphaRhoNew, T);
			double AlphaGNew = AandPAlphaNew[0] + AandPAlphaNew[1]/AlphaRhoNew; 
			
			BetaRhoNew = findDensityBeta(T);
			double[] AandPBetaNew = AandPBeta(BetaRhoNew, T);
			double BetaGNew = AandPBetaNew[0] + AandPBetaNew[1]/BetaRhoNew; 
			
//			System.out.println("********  " + i + " ***********" );
//			System.out.println("Temperature: " + T);
//			System.out.println(AandPAlphaNew[0] + " " + AandPAlphaNew[1] + " "+ AlphaRhoNew+ " "+AlphaGNew);
//			System.out.println(AandPBetaNew[0] + " " + AandPBetaNew[1] + " "+ BetaRhoNew+ " "+BetaGNew);
			
			
//			System.exit(1);
			
//			System.out.println(AandPAlphaNew[0]+" "+AandPAlphaNew[1]);
//			System.out.println(AandPBetaNew[0]+" "+AandPBetaNew[1]);
			deltaGNew = BetaGNew - AlphaGNew;

//			System.exit(1);
			if(deltaGNew > 0.0){
				deltaG0 = deltaGNew;
				T0 = T;
			} else {
				deltaG1 = deltaGNew;
				T1 = T;
			}	
			
			if(i==2000){
				throw new RuntimeException("CANNOT find transition properties!!!!!!");
			}
			++i;
			
		}
		
		this.rhoAlphaTrans = AlphaRhoNew;
		this.rhoBetaTrans = BetaRhoNew;
		
//		System.out.println("rho_alpha: "+AlphaRhoNew + " ;rho_Beta: " + BetaRhoNew);
		return T;
	}
	
	double calcAlphaG(double T){
		double AlphaRho1 = findDensityAlpha(T);
		double[] AandPAlpha1 = AandPAlpha(AlphaRho1, T);
		
		return AandPAlpha1[0] + AandPAlpha1[1]/AlphaRho1; 
		
	}
	
	double calcBetaG(double T){
		double BetaRho1 = findDensityBeta(T);
		double[] AandPBeta1 = AandPBeta(BetaRho1, T);
		
		return AandPBeta1[0] + AandPBeta1[1]/BetaRho1; 
		
	}
	
	int checkNumDensityToUse(double T){
		if(T>49){
			throw new RuntimeException("Temperature for alpha-phase is greater than 49.0!!");
		}
		
		if(T <= 48.0 && T>47.0 ){
			return 4;
		} else if(T <= 47.0 && T>46.0 ){
			return 6;
		} else if(T <= 46.0 && T>45.0 ){
			return 7;
		} else if(T <= 45.0 && T>44.0 ){
			return 9;
		} else if(T <= 44.0 && T>43.0 ){
			return 11;
		} else if(T <= 43.0 && T>42.0 ){
			return 13;
		} else if(T <= 42.0 && T>41.0 ){
			return 15;
		} else {
			return 18;
		}
	}
	
	public static void main (String[] args){

		String fname = "/usr/users/taitan/Research/NormalModes/N2/tp/propagateUncertainty/";
//		String fname = "/usr/users/taitan/Research/NormalModes/N2/tp/findTransT/";
		double tolerance = 1e-6;
		double pGPa = 1e-4;
		
		double [] density = new double[]{0.0229,0.0228,0.0227,0.0226,0.0225,0.0224,0.0223,0.0221,0.022,
				 						 0.0219,0.0218,0.0217,0.0216,0.0215,0.0214,0.0213,0.0212,0.0211};

    	if(args.length > 0){
			pGPa = Double.parseDouble(args[0]);
		}
		double desiredPressure = Pascal.UNIT.toSim(pGPa*1e9);
		FindTransitionTemperature findTransTemp = new FindTransitionTemperature(density, desiredPressure, tolerance, fname);
		
//		for (double t=35.0; t<=46.1; t+=0.1){
//			System.out.println(t + " " + findTransTemp.calcAlphaG(t));
////			System.out.println(t + " " + findTransTemp.calcBetaG(t));
//		}
//		
//		System.exit(1);
		
//		findTransTemp.findDensityBeta(38);
//		System.exit(1);
		
		double TransTemp = findTransTemp.findTransTemperature(38);
		double TransRhoAlpha = findTransTemp.rhoAlphaTrans;
		double TransRhoBeta = findTransTemp.rhoBetaTrans;
		System.out.println("\n********************** Transition Results *********************");
		System.out.println("[I] "+pGPa + " GPa; Transition Temperature found(K): " + TransTemp 
				+ " rhoAlpha: " + TransRhoAlpha + " rhoBeta: " + TransRhoBeta);
		
//		System.exit(1);
		
		 double uAlpha = findTransTemp.UAlpha(TransRhoAlpha, TransTemp);
		 double uBeta = findTransTemp.UBeta(TransRhoBeta, TransTemp);
		 
		 double[] aAlpha= findTransTemp.AandPAlpha(TransRhoAlpha, TransTemp);
		 double[] aBeta= findTransTemp.AandPBeta(TransRhoBeta, TransTemp);
		 
		 double vmAlpha = Constants.AVOGADRO*1e-24/TransRhoAlpha;
		 double vmBeta = Constants.AVOGADRO*1e-24/TransRhoBeta;

		 double TS = (uBeta - uAlpha) - (aBeta[0]-aAlpha[0]);
		 System.out.println("\n[V] vmAlpha (cm^3/mol): "+vmAlpha + " ; vmBeta: " + vmBeta + " ; " + (vmBeta-vmAlpha));
	
		 System.out.println("Tdelta_S: " + TS);
		 System.out.println("Tdelta_S (cal/mol): " + TS*2.390057243);
		 System.out.println("[S] delta_S (cal/mol.K): " + TS*2.390057243/TransTemp);
		 System.out.println("[dS] "+pGPa+" " + TransTemp+" delta_S/delta_V (10^6 kg/(sm^2.K)): " + (TS*2.390057243/TransTemp)/(vmBeta-vmAlpha)*4.184);
		 
	}
	
	protected AkimaSpline akimaSpline;
	protected double[] density; 
	protected double desiredPressure;
	protected double tolerance;
	protected double rhoAlphaTrans, rhoBetaTrans;
	protected double[][] alphaCoeffs, betaCoeffs;
}
