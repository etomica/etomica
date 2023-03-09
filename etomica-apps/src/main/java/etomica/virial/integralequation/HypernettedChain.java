/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.integralequation;

import etomica.math.numerical.SineTransform;


/**
 * 
 * Calculates the hypernetted-chain (HNC) virial coefficients of second to Mth order for any spherically-symmetric Mayer function, fr.
 * 
 * These values differ above third order for the compressilibility and virial routes.
 * 
 * This class has only been tested for the hard sphere and Lennard-Jones potentials.
 * 
 * @author kate
 *
 */

public class HypernettedChain {

	public HypernettedChain() {
	}
	
	public void setRoute(boolean compressibility) {
		this.compressibility=compressibility;
	}
	
	public double[] computeB(double[] fr, int M, int N, double del_r, boolean DCF) {
		
		/*******************************************
		/*******************************************
		 * 
		 * Computes Percus-Yevick approximations of the virial coefficients up to Mth order.  
		 * The second and third coefficients are fully accurate (the PY approximation is exact).
		 * 
		/*******************************************
		/********************************************/
	
		double[] B = new double[M];
		
		double[] fk = SineTransform.forward(fr, del_r);
		
		// Arrays to store the density expansion coefficients of c(r), h(r), c(k) and h(k)
		double[] cnr = new double[N];
		double[] hnr = new double[N];
        double[][] cnk = new double[M][N];
        double[][] hnk = new double[M][N];
		double[][] tnr = new double[M][0];
		
		// Fill zeroth-order density expansion coefficients of c(r) and h(r)
		for (int i=0;i<N;i++) {
			cnr[i] = fr[i];
			hnr[i] = fr[i];
		}
		
		double B2 = -1.0/(2.0)*(fk[0]);
		
		B[0] = B2;
		
		// System.out.println("B2 = " + (B2));
		
		//Compute B3 (from c1) up to BM (from c(M-2))
		for (int m = 1; m <= M-2; m++) {
			
			
			
			/**************************************************************************************
			/**************************************************************************************
			 * Apply the Ornstein-Zernike relation to compute mth-order density expansion of t.
			/**************************************************************************************
			/***************************************************************************************/
			
            cnk[m-1] = SineTransform.forward(cnr, del_r);
			hnk[m-1] = SineTransform.forward(hnr, del_r);
			tnr[m] = OrnsteinZernike.tCompute(cnk, hnk, N, m, del_r);
			
			/**************************************************************************************
			/**************************************************************************************
			 * Apply the Percus-Yevick approximation to compute mth-order density expansion of c.
			/**************************************************************************************
			/***************************************************************************************/
			
			for (int i = 0; i < N; i++) {
				
				
				/*
				 *    m       i
				 *    
				 *    1       1
				 *    
				 *    2       2;  1+1
				 *    
				 *    3       3;  2+1          ; 1+1+1
				 *    
				 *    4       4;  3+1, 2+2     ; 2+1+1;  1+1+1+1
				 *    
				 *    5       5;  4+1, 3+2     ; 3+1+1, 2+2+1; 2+1+1+1;  1+1+1+1+1  
				 *    
				 *    6       6;  5+1, 4+2, 3+3; 4+1+1, 3+
				 */
				
				boolean specialCased = false;
				if (specialCased) {
				
					if (m == 1) {
						hnr[i] = (fr[i]+1) * tnr[1][i];
					} else if (m == 2) {
						hnr[i] = (fr[i]+1) * (tnr[2][i] + (1.0/2.0)*tnr[1][i]*tnr[1][i]);
					} else if (m == 3) {
						hnr[i] = (fr[i]+1) * (tnr[3][i]+tnr[1][i]*tnr[2][i]+(1.0/6.0)*tnr[1][i]*tnr[1][i]*tnr[1][i]);
					} 	
					else if (m == 4) {
						hnr[i] = (fr[i]+1) * (tnr[4][i] + tnr[1][i]*tnr[3][i] + (1.0/2.0)*tnr[2][i]*tnr[2][i] + (1.0/2.0)*tnr[2][i]*tnr[1][i]*tnr[1][i] + (1.0/24.0)*tnr[1][i]*tnr[1][i]*tnr[1][i]*tnr[1][i]);
					} 
				} else {
				
				
					int[] parts = new int[m];
					for (int p=0; p<m; p++) {			
						parts[p]=1;		
					}
					
					//How many ways can we condense m 1's to 1 m?
				    int last = m-1;
				    boolean incomplete = true;
				    hnr[i]=0;
					while (incomplete) {
						
						// count like components that are not zero
						int[] count = new int[m];
						int[] duplicate = new int[m];
						int total = 0;
						for (int p=0; p<m; p++) {  
							if (duplicate[p] == 0) {
							for (int a=p; a<m; a++) {  
								if (parts[a] == parts[p]) {
									count[p] = count[p] + 1;
									total=total+1;
									duplicate[a]=1;
									if (total == m) {
										break;
									}
								} 
							}
							}
							if (total == m) {
								break;
							}
						}
						
						if (i==1){
							
							for (int b=0;b<m; b++) {
								
								System.out.print(parts[b]+"\t");
								
							}
							
							System.out.println();
						}
						
						//add to hnr[m][i]
						
						double plus=1.0;
						for (int p=0; p<m; p++) {  
							if (count[p] != 0 && parts[p] != 0) {	
								double tnra=1.0;
								double afactorial=1.0;
								for (int a=1; a<=count[p]; a++) {
									tnra = tnra*tnr[parts[p]][i];
									afactorial = afactorial*a;
								}
							
								plus = plus *tnra/afactorial;
								
								if (i==1){
									System.out.print(afactorial+"\t");
								}
							}
						}
						hnr[i] += (fr[i]+1)*plus;
						if (i==1){
							System.out.println();
						}
	
						if (parts[0]==m) {
							incomplete = false;
						}
						
						// if incomplete, move 1 from last nonzero element to the leftmost, smallest nonzero element that is not the last nonzero element
						
						int least = m;
						int toID = 0;;
						for (int a=0; a<last; a++) {
							if (parts[a] < least) {
							
								least = parts[a];
								toID = a;
							}
							
							
						}
						
						parts[toID] = parts[toID] + 1;
						parts[last]= parts[last] - 1;
						if (parts[last] == 0) {
							
							last = last -1;
							
						}

					}

					
				}

				cnr[i] = hnr[i] - tnr[m][i];
				
 			}
			
			/*******************************************
			/*******************************************
			 * Calculate (m+2)th virial coefficient
			/*******************************************
			/********************************************/
			
			double[] cmk = SineTransform.forward(cnr, del_r);
			
			double Bm = 0;
			if (compressibility) { //virial route
				Bm = -1.0/(m+2.0)*(cmk[0]); // B3 for m = 1
			} else { //virial route
				
				for (int i=0;i<N-1;i++) {
					double r = del_r*i;
					Bm = Bm + hnr[i]*(rdfdr[i])*r*r;
				}
				
				double r = del_r*(N-1);
				Bm = Bm + 0.5*hnr[N-1]*(rdfdr[N-1])*r*r;
				Bm = Bm*4.0*Math.PI/6.0*del_r;
			}
			
			
			B[m] = Bm;

            //System.out.println("B"+(m+2) + " = "+ (Bm) );
		
			
		
		}
		
		if (DCF) {
			return cnr;
		}
		return B;
	}
	
	public void setrdfdr(double[] rdfdr) {
		this.rdfdr = rdfdr; 
	}
	


	public boolean compressibility = true;
	public double[] rdfdr;
	
}
