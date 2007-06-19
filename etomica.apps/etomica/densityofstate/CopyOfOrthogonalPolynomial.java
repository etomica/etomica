package etomica.densityofstate;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;


public class CopyOfOrthogonalPolynomial {
	public static void main (String[] args){
	
		double[] t = new double [args.length];
		double[][] p = new double [args.length] [0];
		double[][] u = new double [args.length] [0];
		double[][] b = new double [args.length] [7];
		double[][] lnDOS = new double [args.length] [0];
		
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
	            lnDOS[i] = new double[numLines-1];
	           
	            t[i] = Double.parseDouble(bufReader.readLine());
	            //System.out.println(t[i]);
	            
	            for (int j = 0; j < numLines-1; j++){
	            	String line = bufReader.readLine();
	            	String [] energyProb = line.split("[\t ]+");
	            	u[i][j] = Double.parseDouble(energyProb[0]);
	            	p[i][j] = Double.parseDouble(energyProb[1]);
	            	
	            	//System.out.println(u[i][j] + " " + p[i][j]);
	            }
	            
	        }catch (IOException e){
	        	e.printStackTrace();
	            throw new RuntimeException("Cannot read "+args[i]+", caught IOException: " + e.getMessage());
	        }
		
	        for (int n = 0; n < 7; n++ ){                     //different n of phi
	        	double numer = 0;
	        	double denom = 0;
	        	
	        	for (int j = 0; j < u[i].length; j++){         //number of lines of data 
	        		double phinj = phi(n,u[i][j]);
	        		numer += (Math.log(p[i][j])-(u[i][j]/t[i]))*phinj;
	        		denom += phinj*phinj;
	        	}
	        	
	        	b[i][n] = numer/denom;							//calculate for bn
	        	System.out.println(b[i][n]);
	        }
		
	        System.out.println("\n");
	        //calculate for DOS at each energy
	        
	        for(int j = 0; j < u[i].length; j++){
	        	for (int n = 0; n < 7; n++){
	        		lnDOS[i][j] += b[i][n]*phi(n,u[i][j]);
	        		
	        		//System.out.print(phi(n,t[i],u[i][j])+ " ");
	        		}
	        	//System.out.println(" ");
	        	System.out.println(u[i][j] + " " + lnDOS[i][j] + " " + Math.exp(lnDOS[i][j]) + " " + Math.exp(lnDOS[i][j])*Math.exp(-u[i][j]/t[i]));
	        	
	        }
	    	
		}
	}
		
	public static double phi (int i, double U){
		double U2 = U*U;
		double U4 = U2*U2;
		
		switch (i){
		case 0: return 1;
		case 1: return U;
		case 2: return (-1/2) + (3/2)*U2;
		case 3: return (-3/2) + (5/2)*U2*U;
		case 4: return (3/8) - (30/8)*U2 + (35/8)*U4;
		case 5: return (15/8)*U -(70/8)*U2*U + (63/8)*U4*U;
		case 6: return (-5/16) +(105/16)*U2 - (315/16)*U4 + (213/16)*U4*U2;
		
		default: throw new IllegalArgumentException("i must be less than 7");
		} 
	}
	
	
}
	