//Input Temperature has to be in increasing order (top to bottom)


package etomica.densityofstate;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;


public class FerrenbergSwendsen {
	public static void main (String[] args){
	
		double[] t = new double [args.length];
		double[][] p = new double [args.length] [0];
		double[][] u = new double [args.length] [0];
		double[] z = new double [args.length];
		double deltaU = Double.NaN;
		double uMax = Double.NEGATIVE_INFINITY;
		double uMin = Double.POSITIVE_INFINITY;
		
		for (int i = 0; i < args.length; i++){
			FileReader fileReader;
	        try {
	            fileReader = new FileReader(args[i]);
	        }catch(IOException e) {
	            throw new RuntimeException("Cannot open "+args[i]+", caught IOException: " + e.getMessage());
	        }
	        try {
	            BufferedReader bufReader = new BufferedReader(fileReader);
	            int numLines=0;
	            while (bufReader.readLine()!=null){
	            	numLines++;
	            }	           
	            bufReader.close();
	            fileReader.close();
	            fileReader = new FileReader (args[i]);
	            bufReader = new BufferedReader(fileReader);
	            p[i] = new double[numLines-1];
	            u[i] = new double[numLines-1];
	           
	            t[i] = Double.parseDouble(bufReader.readLine());
	            //System.out.println(t[i]);
	            
	            for (int j = 0; j < numLines-1; j++){
	            	String line = bufReader.readLine();
	            	String [] energyProb = line.split("[\t ]+");
	            	u[i][j] = Double.parseDouble(energyProb[0]);
	            	p[i][j] = Double.parseDouble(energyProb[1]);
	            	//System.out.println(u[i][j] + " " + p[i][j]);
	            }
	            
	            if (i == 0) {
	            	deltaU = u[i][1] - u[i][0]; 	//determine bin size
	            }
	            
	            if (uMin > u[i][0]){
	            	uMin = u[i][0];					//determine minimum energy level
	            }
	            
	            if (uMax < u[i][numLines-2]){
	            	uMax = u[i][numLines-2];		//determine maximum energy level
	            }
	            
	        }catch (IOException e){
	        	e.printStackTrace();
	            throw new RuntimeException("Cannot read "+args[i]+", caught IOException: " + e.getMessage());
	        }
		}
		
		int numBin = (int) Math.round(((uMax - uMin) / deltaU) + 1);
		int[] binOffset = new int [p.length]; 
	
		for (int k = 0; k < p.length; k++){
			binOffset[k] = (int) Math.round((u[k][0] - uMin) / deltaU);			//to determine the bin difference of the   
			z[k] = 1;															//first bin of each histograms relative to uMin
		}
		
		double[] total = new double [p.length];
		boolean converge = false;
		
		while (!converge){
			
			for (int j = 0; j < numBin  ; j++ ){
				
				double pSum = 0;
				double denomSum = 0;
				
				double sameEnergy = uMin + j * deltaU;							//each energy level ranging uMin to uMax
				
				for (int i = 0; i < p.length; i++){
					
					int index = j - binOffset[i];
					denomSum += Math.exp(-(1/ t[i]) * sameEnergy) /z[i];
					
					
					if (index < 0 || index > p[i].length - 1){
						continue;	
					}
					
					pSum += p[i][index];
				}
							
				for (int i = 0; i < p.length; i++){
					
					total[i] += deltaU * Math.exp(-(1/ t[i]) * sameEnergy) * pSum / denomSum;
					
					//if (Double.isNaN(total[i])){
					//	throw new RuntimeException(" Oh no!!");
					//}
				
				}
			}
			
			converge = true;
			
			for (int i = 0; i < p.length ; i++){
				total[i] = total[i] / total[p.length-1];
				
				if (Math.abs((total[i] - z[i])/ z[i]) > 0.00001){
					 converge = false;
				}
				
				z[i] = total[i];
				System.out.print(z[i] + " ");		//larger z[i] is set to 1
				total[i] = 0;							//In decending order
			}
			System.out.println(" ");	
		}
	
		//Combine all into master histogram
		
		for (int j = 0; j < numBin  ; j++ ){
			
			double pSum = 0;
			double denomSum = 0;
			
			double sameEnergy = uMin + j * deltaU;							//each energy level ranging uMin to uMax
			
			for (int i = 0; i < p.length; i++){
				
				int index = j - binOffset[i];
				denomSum += Math.exp(-(1/ t[i]) * sameEnergy) /z[i];
				
				if (index < 0 || index > p[i].length - 1){
					continue;	
				}
				
				pSum += p[i][index];
			}
						
			System.out.println(sameEnergy + " " + pSum/denomSum);
			
		}
		
	
	}
}
	