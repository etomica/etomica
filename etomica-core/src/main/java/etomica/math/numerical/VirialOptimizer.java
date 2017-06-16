/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.numerical;

import etomica.util.ParameterBase;
import etomica.util.ReadParameters;

/**
 * 
 * Class to determine B9 value that best fit the simulation results after 
 *   Pade Approximation [K/L]  
 * 
 * @author Tai Boon Tan
 *
 */
public class VirialOptimizer {
	
	public VirialOptimizer(String filenameP, String filenameB){
		
		pSimData = ArrayReader1D.getFromFile(filenameP);
		bVirialFromFile = ArrayReader1D.getFromFile(filenameB);
		bVirial = new double[bVirialFromFile.length + 1];
		rho = new double[pSimData.length];
		pPade = new double[pSimData.length];
		
		for (int i=0; i<bVirialFromFile.length;i++){
			bVirial[i] = bVirialFromFile[i][0];
		}
		
		for (int i=0; i<pSimData.length;i++){
			rho[i] = pSimData[i][0];
		}
		
	}
	
	public void calcpPade(int K, int L, double bGuess){
		/*
		 * Pade Approximation
		 */
		double pPadeNum, pPadeDenom;
		
		/*
		 * Optimizing higher-order B-Coefficient 
		 */
		bVirial[bVirial.length-1] = bGuess;
		
		PadeApproximation pade = new PadeApproximation(bVirial, K, L);
		pade.solveCoefficients();
		double[] A = pade.getA();
		double[] C = pade.getB();
		
		/*
		 * Pressure Calculation from Pade
		 */
		
		for (int irho=0; irho< rho.length; irho++){
			pPadeNum = 0.0;
			pPadeDenom = 0.0;
						
			for (int i=0; i<A.length;i++){
				pPadeNum += A[i]*Math.pow(rho[irho],i);
			}
			
			for (int i=0; i<C.length;i++){
				pPadeDenom += C[i]*Math.pow(rho[irho],i);
			}
			
			pPade[irho] = pPadeNum/pPadeDenom;
			
		}
	}
	
	public double minimumScreening(int K, int L, double min, double max){
		
		int nPoint = 50000;
		double globalMin = 50.0; //this is arbitrary value
		double interval = (max-min)/nPoint;
		double b = min;
		double bOpt = min;
		
		for (int i=0; i<nPoint; i++){
			
			calcpPade(K, L, b);
			
			double TSS = 0.0;
			for (int j=0; j<pPade.length; j++){
				diff = (pPade[j]-pSimData[j][1]);
				TSS += diff*diff;
			}
			if(TSS < globalMin){
				globalMin = TSS;
				bOpt = b;
			}
			b += interval;

		}
		
		return bOpt;
	}
	
	public void optimizeHigherbVirial(int K, int L, double minb, double maxb){
		
		int bootstrap = 0;
		double[] allTSS = new double[3];
		double[] allb = new double[3];
		
		double bGuess = minb;
		int counter = 0;
		
		while (true){
			
			
			calcpPade(K, L, bGuess);
			
			double TSS = 0.0;
			for (int i=0; i<pPade.length; i++){
				diff = (pPade[i]-pSimData[i][1]);
				TSS += diff*diff;
				
			}
			
			if(bootstrap < 3){
				allb[bootstrap] = bGuess;
				allTSS[bootstrap] = TSS;
				bootstrap++;
				bGuess += 0.5*(maxb-minb);
			}
		    else {
                if (bGuess > allb[2]) {
                    allb[0] = allb[1];
                    allb[1] = allb[2];
                    allb[2] = bGuess;
                    allTSS[0] = allTSS[1];
                    allTSS[1] = allTSS[2];
                    allTSS[2] = TSS;
                }
                else if (bGuess < allb[0]) {
                    allb[2] = allb[1];
                    allb[1] = allb[0];
                    allb[0] = bGuess;
                    allTSS[2] = allTSS[1];
                    allTSS[1] = allTSS[0];
                    allTSS[0] = TSS;
                }
                else if (allTSS[2] > allTSS[0]) {
                    maxb = allb[2];
                    if (bGuess > allb[1]) {
                        allTSS[2] = TSS;
                        allb[2] = bGuess;
                    }
                    else {
                        allTSS[2] = allTSS[1];
                        allb[2] = allb[1];
                        allTSS[1] = TSS;
                        allb[1] = bGuess;
                    }
                }
                else {
                    minb = allb[0];
                    if (bGuess < allb[1]) {
                        allTSS[0] = TSS;
                        allb[0] = bGuess;
                    }
                    else {
                        allTSS[0] = allTSS[1];
                        allb[0] = allb[1];
                        allTSS[1] = TSS;
                        allb[1] = bGuess;
                    }
                }

            }
		
	         if (bootstrap == 3) {
	        	 
	                double dc01 = allb[1]-allb[0];
	                double dc12 = allb[2]-allb[1];
	                double du01 = allTSS[1]-allTSS[0];
	                double du12 = allTSS[2]-allTSS[1];
	                double dudc01 = du01/dc01;
	                double dudc12 = du12/dc12;
	                double m = (dudc12-dudc01)/(0.5*(dc01+dc12));
	                bGuess = 0.9*(0.5*(allb[1]+allb[2]) - dudc12/m) + 0.1*(0.5*(allb[0]+allb[2]));
	                if (bGuess == allb[1] || bGuess == allb[2]) {
	                    bGuess = 0.5*(allb[1] + allb[2]);
	                }
	                if (bGuess == allb[0] || bGuess == allb[1]) {
	                    bGuess = 0.5*(allb[1] + allb[0]);
	                }
	                if (bGuess < minb) {
	                    bGuess = 0.5*(minb + allb[0]);
	                }
	                if (bGuess > maxb) {
	                    bGuess = 0.5*(maxb + allb[2]);
	                }
	                        
	                if (bGuess == allb[0] || bGuess == allb[1] || bGuess == allb[2]) {
	                   	break;
	                }
	                ++counter;
	                
	                if(counter > 1e5){
	                	System.out.println("<Java.VirialOptimizer> Could not find minimum!!");
	                	break;
	                }
	            }
			
		}
				
		System.out.println(bGuess);
	}
	
	
	public static void main(String[] args){
		
		VirialParam params = new VirialParam();
	    String inputFilename = null;
        if (args.length > 0) {
            inputFilename = args[0];
        }
        if (inputFilename != null) {
            ReadParameters readParameters = new ReadParameters(inputFilename, params);
            readParameters.readParameters();
        }
        
		String filenameP = params.filenameP;
		String filenameB = params.filenameB;
		int K = params.K;
		int L = params.L;
		
		VirialOptimizer vOpt = new VirialOptimizer(filenameP, filenameB);
			
		double x = vOpt.minimumScreening(K, L, -5, 5);
		double min = (1-0.05)*x;
		double max = (1+0.05)*x;
		
		vOpt.optimizeHigherbVirial(K, L, min, max);		
	
	}
	
	protected double diff;
	protected double[] bVirial, rho, pPade;
	protected double[][] pSimData, bVirialFromFile;
	protected double bGuess;
	
	public static class VirialParam extends ParameterBase {
	        public String filenameP = "/tmp/foo";
	        public String filenameB = "/tmp/Bn9"; 
	        public int K = 5;
	        public int L = 3;
	}
}
