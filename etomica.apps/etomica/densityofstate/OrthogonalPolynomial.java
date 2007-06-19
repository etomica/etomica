package etomica.densityofstate;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;


public class OrthogonalPolynomial {
	public static void main (String[] args){
	
		double[] t = new double [args.length];
		double[][] p = new double [args.length] [0];
		double[][] u = new double [args.length] [0];
		double[][] b = new double [args.length] [10];
		double[][] DOS = new double [args.length] [0];
		
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
	            DOS[i] = new double[numLines-1];
	           
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
		
	        for (int n = 0; n < 10; n++ ){                     //different n of phi
	        	double numer = 0;
	        	double denom = 0;
	        	
	        	for (int j = 0; j < u[i].length; j++){         //number of lines of data 
	        		
	        		double phinj = phi(n, t[i],u[i][j]);
	        		double newphinj = 1;
	        		
	        		for (int k = 0; k < 30; k++){
	        			newphinj *= phinj;
	        			k++;
	        		}
	        		
	        		numer += p[i][j]*newphinj;
	        		denom += newphinj*newphinj*Math.exp(-u[i][j]/t[i]);
	        	}
	        	
	        	b[i][n] = numer/denom;							//calculate for bn
	        	System.out.println(b[i][n]);
	        }
		
	        System.out.println("\n");
	        //calculate for DOS at each energy
	        
	        for(int j = 0; j < u[i].length; j++){
	        	for (int n = 0; n < 10; n++){
	        		DOS[i][j] += b[i][n]*phi(n,t[i],u[i][j]);
	        		
	        		//System.out.print(phi(n,t[i],u[i][j])+ " ");
	        		}
	        	//System.out.println(" ");
	        	System.out.println(u[i][j] + " " + DOS[i][j] + " " + DOS[i][j]*Math.exp(-u[i][j]/t[i]));
	        	
	        }
	    	
		}
	}
		
	public static double phi (int i, double T, double U){
		double UT = U/T;
		double UT2 = UT*UT;
		double UT4 = UT2*UT2;
		double UT6 = UT4*UT2;
		double UT8 = UT4*UT4;
		
		switch (i){
		case 0: return 1/Math.sqrt(T);
		case 1: return (1-  UT)/Math.sqrt(T);
		case 2: return (1-2*UT+ (1/2)*UT2)/Math.sqrt(T);
		case 3: return (1-3*UT+ (3/2)*UT2- (1/6)*UT*UT2)/Math.sqrt(T);
		case 4: return (1-4*UT+   3  *UT2- (2/3)*UT*UT2+ (1/24)*UT4)/Math.sqrt(T);
		case 5: return (1-5*UT+   5  *UT2- (5/3)*UT*UT2+ (5/24)*UT4-(1/120)*UT*UT4)/Math.sqrt(T);
		case 6: return (1-6*UT+(15/2)*UT2-(10/3)*UT*UT2+ (5/8) *UT4- (1/20)*UT*UT4+(1/720)*UT6)/Math.sqrt(T);
		case 7: return (1-7*UT+(21/2)*UT2-(35/6)*UT*UT2+(35/24)*UT4- (7/40)*UT*UT4+(7/720)*UT6-(1/5040)*UT*UT6)/Math.sqrt(T);
		case 8: return (1-8*UT+  14  *UT2-(28/3)*UT*UT2+(35/12)*UT4- (7/15)*UT*UT4+(7/180)*UT6- (1/630)*UT*UT6+(1/40320)*UT8)/Math.sqrt(T);
		case 9: return (1-9*UT+  18  *UT2-  14  *UT*UT2+ (21/4)*UT4-(21/20)*UT*UT4+( 7/60)*UT6- (1/140)*UT*UT6+ (1/4480)*UT8-(1/362880)*UT*UT8)/Math.sqrt(T);
		
		default: throw new IllegalArgumentException("i must be less than 10");
		} 
	}
	
	
}
	