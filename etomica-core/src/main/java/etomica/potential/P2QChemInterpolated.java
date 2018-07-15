/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.util.ParameterBase;
import etomica.math.numerical.AkimaSpline;

/**
 * Pair potential for argon interpolated from Q-Chem results.  
 * 
 * In this class, only the pair potential is valid, not the gradients, etc.  I am unlikely to ever include those...
 * 
 * To use this class, one must call initialize(), setDisp(), setSCF()
 * 
 * If one sets disp to be true with setDisp(), one must also call setDampingParamters();
 *
 * @author Kate Shaul
 */
public class P2QChemInterpolated extends Potential2SoftSpherical {
    
    public P2QChemInterpolated(Space space) {
    	
        super(space);
   
    }

    /**
     * The energy u.
     */
    public void setDampingParams(int a1, int a2, double Rvdw, int basis, boolean fixedRvdw) {
    	this.a1 =a1;
    	this.a2 = a2;
    	this.Rvdw = Rvdw;
    	this.basis = basis;
    	this.fixedRvdw = fixedRvdw;
    }
    
    public void setDisp(boolean disp) {
    	this.disp = disp;
    }
    
    public void setSCF(boolean scf) {
    	this.scf = scf;
    }
   
    public void initialize() {
    	getData();
    	getTheta();
    	getSplines();	
    }
    
    public double u(double r2) {
    	
    	double r = Math.sqrt(r2);
    	
    	double uDisp = 0;
    	double uSCF = 0;
    	
    	if (r < 1.6) {
    		return Double.POSITIVE_INFINITY;
    	}
    	
    	if (disp) {
    		uDisp = getUDisp(r);
    		//System.out.println(r+ "  " +uDisp*1e6);
    		if (r > 8) {
    			return uDisp*JPerHartree/k;
    		}
    	} 
    	
    	if (scf) {
    		uSCF = getUSCF(r);
    		//System.out.println(r+ "  " +uSCF*1e6);
    	} 

    	//System.out.println(r+ "  " +(uDisp+uSCF)*1e6);
    	return (uDisp+uSCF)*JPerHartree/k; // K;
    }
   
    public void getSplines() {
		
		double[] ratio = new double[rData.length];
		
		for (int i=0;i<rData.length;i++) {
			
			double f = theta[0]*Math.exp(theta[1]*rData[i]*rData[i]);
			
			// The ratio that will be interpolated
			ratio[i]=uSCFData[i]/f;
			
		}
		
		splineRatio.setInputData(rData, ratio);
		splineUSCF.setInputData(rData, uSCFData);
		splineC6.setInputData(r2Data, C6Data);
		splineC8.setInputData(r2Data, C8Data);
		splineC10.setInputData(r2Data, C10Data);
		splineRc.setInputData(r2Data, RcData);
    	
    	
    }
   
    
    public double getUSCF(double r) {
    	
    	int N = rData.length;
    	double uSCF=0;
    	int i = 0;
    	
    	if (r < 0.1) {
    		return Double.POSITIVE_INFINITY;
    	} else {
    		
    		for (int n=0; n<N; n++) {	
    			if (r == rData[n]) {
    				i = n;
    			}
    		}
    			
    		if (i != 0) {
    			
    			uSCF = uSCFData[i];
    			
    		} else {
    		
				double[] rA = new double[] {r};
				
    			if (r < 3.5) {
    				
    				double[] ratio = splineRatio.doInterpolation(rA);
    				
    				double f = theta[0]*Math.exp(theta[1]*r*r);
    				
    				uSCF = f*ratio[0];
	   			
    			} else {
    				
    	    		double[] energyA = splineUSCF.doInterpolation(rA);
    	    		uSCF = energyA[0];
    	    		
    			}
			}

    	}
    	
		return uSCF;
    }
    
    public double getUDisp(double r) {
    	
    	//////////////////////////////////////////////////////////////////////////////////
    	// supplemental dispersion component of Becke-Johnson model
    	//////////////////////////////////////////////////////////////////////////////////
    	
    	
    	double uDisp=0;
    	double C6=0;
		double C8=0;
		double C10=0;
		double Rc=0;
		int N = r2Data.length;
    	int i = 0;
		
    	if (r < 0.1) {
    		return Double.POSITIVE_INFINITY;
    	} else {
 
    		for (int n=0;n<N; n++) {
    			
    			if (r == r2Data[n]) {
    				i = n;
    			}	
    			
    			/*if (r2Data[i]==10) {
    				System.out.println(i);
    				System.exit(0);
    			}*/
    			
    		}
    		
    		if (r > 10) {
    			
    			//n = 102; // 8 Angstroms
    			i = 122; // 10 Angstroms
    			//n = 124; // 12 Angstroms
	    		
    		} 
    		
    		if (i != 0) {

				C6 = C6Data[i];
				C8 = C8Data[i];
				C10 = C10Data[i];
				Rc = RcData[i];
    			
    		} else {
	    		
    			double[] rA = new double[] {r};

	    		double[] C6A = splineC6.doInterpolation(rA);
	    		C6 = C6A[0];  		

	    		double[] C8A = splineC8.doInterpolation(rA);
	    		C8 = C8A[0]; 		

	    		double[] C10A = splineC10.doInterpolation(rA);
	    		C10 = C10A[0];

	    		double[] RcA = splineRc.doInterpolation(rA);
	    		Rc = RcA[0];

    		}
    	
    			
    		if (!fixedRvdw) {
        		Rvdw = (a1*Rc + a2)/100;
        	}
			double energy6  =  -C6/(Math.pow(Rvdw,6) + Math.pow(r,6)) *Math.pow(0.529177209,6);
			double energy8  =  -C8/(Math.pow(Rvdw,8) + Math.pow(r,8)) *Math.pow(0.529177209,8);
			double energy10 = -C10/(Math.pow(Rvdw,10)+ Math.pow(r,10))*Math.pow(0.529177209,10);
			uDisp = energy6 + energy8 + energy10; //Hartrees
			
			return uDisp;

    	}

    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        
        return 0;
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
     
        return 0;
    }
            
    /**
     *  Integral used for corrections to potential truncation.
     */
    public double uInt(double rC) {
        
        return 0;  //complete LRC is obtained by multiplying by N1*N2/V
    }
    
    

    public static void main(String[] args)  {
    	
    	Space space = Space3D.getInstance();
    	P2QChemInterpolated p2 = new P2QChemInterpolated(space);
    	
		DampingParams params = new DampingParams();
	     
		if (args.length == 5 ) {
			params.a1 = Integer.parseInt(args[0]);
			params.a2 = Integer.parseInt(args[1]);
			params.RvdwF = Double.parseDouble(args[2]);
			params.basis = Integer.parseInt(args[3]);
			params.fixedRvdw = Boolean.parseBoolean(args[4]);
	    } 
		
		int a1 = params.a1;
		int a2 = params.a2;
		double RvdwF = params.RvdwF;
		int basis = params.basis;
		boolean fixedRvdw = params.fixedRvdw;
		p2.setDampingParams(a1,a2,RvdwF, basis, fixedRvdw);
		p2.setDisp(true);
		p2.setSCF(true);
		p2.initialize();
		
		
		//double r = 2.28254;
		double r = 3.757178;
		System.out.println(r+ "  " +p2.u(r*r)*k/JPerHartree*1e6);
		/*
		for (int i=200;i<=500;i++) {
			double r = 0.01*i;
			double u=p2.u(r*r);
			System.out.println(r+ "  " +u*k/JPerHartree);
			
		} 
		*/
		
    	
    	double[] min = p2.minimum(p2);
    	
    	//double Rvdw = p2.optimizeRvdwWMin(p2);
    	
    	//double RvdwuTT = p2.optimizeRvdwWuTT3757(p2);
    	
    	
    	
    	
    	
    	//r=3.757178;
    	//double r=3.757;
    	//System.out.println("r12=	"+r + ", u12 =    " +p2.u(r*r)*k/JPerHartree*1e6);
    	
    	
    }
   
    public double[] minimum(P2QChemInterpolated p2) {
    	double r = 3.6;
    	double umin=0;
    	double rmin=0;
    	
    	while (r<3.9) {
    		r = r + 0.0001;
    		double u = p2.u(r*r)*k/JPerHartree*1e6; 
    		
    		if (u < umin) {
    			umin = u;
    			rmin =r;
    		}
    		
    		System.out.println(r+"  "+u);
    	}
    	
    	System.out.println(Rvdw+"  "+rmin+"  "+umin);
    	
    	return new double[] {rmin, umin};
    }
    
    public double optimizeRvdwWMin(P2QChemInterpolated p2) {
    	
    	Rvdw=3.55;
        double error=100;
        double uminF = 0;
    	double rminF = 0;
     	while (Rvdw<=3.85) {
     		
     		Rvdw = Rvdw + 0.01;
     		p2.setDampingParams(a1,a2,Rvdw, basis, fixedRvdw);
     		double[] min = p2.minimum(p2);
     		double rmin = min[0];
     		double umin = min[1];
     		double e1 =(umin+454)/454;
	    	double e2 = 3*(rmin-3.757)/3.757;
	    	double e = e1*e1+e2*e2;
	    	//double e = e1*e1;
	    	
	    	if ( e < error) {
	    		RvdwF = Rvdw;
	    		uminF = umin;
	    		rminF = rmin;
	    		error = e;
	    	}
	    	
	    	 
	    	//System.out.println(Rvdw+"	"+rmin+ "    " +umin+ "    " +e);
	    	
    	}
     	
     	System.out.println(Rvdw+"	"+rminF+ "    " +uminF+ "    " +error);
    	
    	return RvdwF;
    }
    
    
    
    public double optimizeRvdwWuTT3757(P2QChemInterpolated p2) {
    	
    	double r=3.757;
    	
    	Rvdw=3.55;
        double error=100;
        double uF = 0;
     	while (Rvdw<=3.85) {
     		
     		Rvdw = Rvdw + 0.01;
     		p2.setDampingParams(a1,a2,Rvdw, basis, fixedRvdw);
     		double u = p2.u(r*r)*k/JPerHartree*1e6;
    		if (Math.abs(u+454) < error) {
    			RvdwF=Rvdw;
    			uF = u;
    			error = Math.abs(u+454);
    		}
    		System.out.println(Rvdw+"  "+r+"  "+u);
	    	
    	}
     	
     	System.out.println(Rvdw+ "    " +uF+ "    " +error);
    	
    	return RvdwF;
    }
    

    
    
	public static class DampingParams extends ParameterBase {
		/*
		protected int a1 = 79;	        
	    protected int a2 = 136;   
	    protected int basis = 3;
	    protected double RvdwF = 3.70;   
	    */
	    protected int a1 = 80;	        
	    protected int a2 = 149;   
	    protected int basis = 2;
	    protected double RvdwF = 3.688; 
	    
	    //protected int basis = 4;
	     
	    protected boolean fixedRvdw = false; 
	   
	    
	 

	}
	
	public void getTheta() {
		
		////////////////////////////////////////////////////////////
		// Remove x and y elements that don't make sense.
		////////////////////////////////////////////////////////////
		
		double[] xA = rData;
		double[] yA = uSCFData;
		
		int N = yA.length;
		int min = 0; 
		int max = N;
		for (int n = 0; n < N; n++) {	
			//if (yA[n] > 1000) {	
			if (xA[n] < 2) {	
				min = n;
			}
			if (xA[n] < 3.5) {	
				max = n+1;
				
			}
			//System.out.println(xA[n]+"  "+yA[n]);
		} 
		
		N = max - min;
		
		//System.out.println(min+" "+max+" "+N);
		
		double [] xB = new double[N]; 
		double [] yB = new double[N];
	
		for (int n = 0; n < N; n++) {	
			xB[n] = xA[n+min];
			yB[n] = yA[n+min];
			
		}
		
		////////////////////////////////////////////////////////////
		// Prepare and solve normal equations for X = [ln(theta0); theta1]
		// A'*A*X = A'*B -->  X = (A'*A)^-1*A'*B
		//
		// y = theta0*exp(theta1*x) 
		// ln(y) = ln(theta0) + theta1*x
		////////////////////////////////////////////////////////////
		
		double [][] a = new double[N][2]; 
		double [] b = new double[N];
		
		for (int n = 0; n < N; n++) {	
			a[n][0] = 1;
			a[n][1] = xB[n]*xB[n];
			b[n] = Math.log(yB[n]);
		}
		
//		Matrix A = new Matrix(a);
//		Matrix B = new Matrix(b, N);
//
//		Matrix At = A.transpose();
//		Matrix X = ((At.times(A)).inverse()).times(At.times(B));
			
		double theta0 = 0;
		double theta1 = 0;
		
		////////////////////////////////////////////////////////////
		// Normal equations
		// A'*A*X = A'*B -->  X = (A'*A)^-1*A'*B
		//
		// y = theta0*exp(theta1*x) 
		// ln(y) = ln(theta0) + theta1*x
		////////////////////////////////////////////////////////////
		
		double[] f = new double[N];
		double[] e = new double[N];
		double[] f2 = new double[N];
		double[] e2 = new double[N];
		
		for (int n=0; n<N; n++) {
						
			f[n] = theta0*Math.exp(theta1*(xB[n]*xB[n]));
			
			f2[n] = Math.log(theta0) + theta1*xB[n]*xB[n];
			
			e[n] = f[n] - yB[n];
			
			e2[n] = f[n] - b[n];

			//System.out.println(xB[n]+"  "+yB[n]+"  "+f[n]);
			
		}
		
//		Matrix E = new Matrix(e, N);
//		double error = E.norm2();
//
//		Matrix E2 = new Matrix(e2, N);
		double error2 = 0;
			
		theta = new double[] {theta0,theta1};
		//System.out.println("theta0 = " + theta0 + ", theta1 = "+ theta1 + ", error = " + error+ ", error2 = " + error2);
			
		//System.exit(0);

	}
	
	 public void getData() {
		 	String d;
	    	if (basis == 4) {
	    		d = "/usr/users/kate/ArData/PW86PBE/augccpvqz/";
	    	} else if (basis == 3) {
	    		d = "/usr/users/kate/ArData/PW86PBE/augccpvtz/";
	    	} else if (basis == 2) {
	    		d = "/usr/users/kate/ArData/PW86PBE/augccpvdz/";
	    	} else {
				throw new RuntimeException("Specify different basis.");
			}
				
	
			double[] dataUSCF = new double[1000];
			double[] dataR = new double[1000];
			
			int count = 0;
			try{			    			
    			FileReader fileReader = new FileReader(d+"u12NoDisp.dat");			
    			BufferedReader bufReader = new BufferedReader(fileReader);
    			String line;
    			
    			while ((line = bufReader.readLine()) != null) {
    				
    				dataUSCF[count] = Double.parseDouble(line.split(" +")[1]); 
    				dataR[count] = Double.parseDouble(line.split(" +")[0]); 
    				count++;
    			}
    		}catch (IOException e){	
    			throw new RuntimeException(e);
    		}
    		
    		uSCFData = new double[count];
    		rData = new double[count];
    		for (int n=0;n<count;n++) {
    			uSCFData[n] = dataUSCF[n];
    			rData[n] = dataR[n];
    			//System.out.println(rData[n]+" "+uSCFData[n]);
    		}
    		//System.exit(0);
			
		
    		// Clean: 
			for (int i = 0; i< uSCFData.length; i++) {
				if (uSCFData[i] > 1.05e+03) {
					uSCFData[i] = Double.POSITIVE_INFINITY;
		    	}
			}
				
			double[] dataR2 = new double[1000];	
			double[] dataC6 = new double[1000];
			count = 0;
			try{			    			
    			FileReader fileReader = new FileReader(d+"C6.dat");			
    			BufferedReader bufReader = new BufferedReader(fileReader);
    			String line;
    			while ((line = bufReader.readLine()) != null) {	
    				dataC6[count] = Double.parseDouble(line.split(" +")[1]);  
    				dataR2[count] = Double.parseDouble(line.split(" +")[0]);
    				count++;
    			}
    		}catch (IOException e){	throw new RuntimeException(e);}
			
    		double[] dataC8 = new double[1000];
			count = 0;
			try{			    			
    			FileReader fileReader = new FileReader(d+"C8.dat");			
    			BufferedReader bufReader = new BufferedReader(fileReader);
    			String line;
    			while ((line = bufReader.readLine()) != null) {	
    				dataC8[count] = Double.parseDouble(line.split(" +")[1]);  
    				count++;
    			}
    		}catch (IOException e){	throw new RuntimeException(e);}
    		
    		double[] dataC10 = new double[1000];
    		count = 0;
			try{			    			
    			FileReader fileReader = new FileReader(d+"C10.dat");			
    			BufferedReader bufReader = new BufferedReader(fileReader);
    			String line;
    			while ((line = bufReader.readLine()) != null) {	
    				dataC10[count] = Double.parseDouble(line.split(" +")[1]);  
    				count++;
    			}
    		}catch (IOException e){	throw new RuntimeException(e);}
    		
    		double[] dataRc = new double[1000];
    		count = 0;
			try{			    			
    			FileReader fileReader = new FileReader(d+"Rc.dat");			
    			BufferedReader bufReader = new BufferedReader(fileReader);
    			String line;
    			while ((line = bufReader.readLine()) != null) {	
    				dataRc[count] = Double.parseDouble(line.split(" +")[1]);  
    				count++;
    			}
    		}catch (IOException e){	throw new RuntimeException(e);}
			
    		C6Data  = new double[count];
    		C8Data  = new double[count];
    		C10Data = new double[count];
    		RcData  = new double[count];
    		r2Data  = new double[count];
    		for (int n=0;n<count;n++) {
    			C6Data[n]  = dataC6[n];
    			C8Data[n]  = dataC8[n];
    			C10Data[n] = dataC10[n];
    			RcData[n]  = dataRc[n];
    			r2Data[n]  = dataR2[n];
    		}
    	
	    }
   
	protected int a1;	        
    protected int a2;   
    protected int basis;
    protected double RvdwF;   
    protected static double Rvdw; 
    protected boolean fixedRvdw ; 
    public boolean disp = true;
    public boolean scf = true;
    protected static double k = 1.3806503e-23; //J/K   1.3806503e-23
    protected static double JPerHartree = 4.359744e-18;  //4.359744e-18
    protected static double Na = 6.0221415e23;  
    protected boolean dense = false;
    public double[] theta;
    protected AkimaSpline splineRatio = new AkimaSpline();
    protected AkimaSpline splineUSCF = new AkimaSpline();
    protected AkimaSpline splineC6 = new AkimaSpline();
    protected AkimaSpline splineC8 = new AkimaSpline();
    protected AkimaSpline splineC10 = new AkimaSpline();
    protected AkimaSpline splineRc = new AkimaSpline();
    protected double[] rData; 
    protected double[] uSCFData; 
    protected double[] C6Data; 
    protected double[] C8Data; 
    protected double[] C10Data; 
    protected double[] RcData; 
    protected double[] r2Data; 
    
    
    }
