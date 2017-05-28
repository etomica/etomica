/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.math.SpecialFunctions;
import etomica.units.Electron;
import etomica.units.Kelvin;

/**
 * two points cluster diagrams calculation
 * Integral I2,S2,mean activity coefficient,osmotic coefficient of ionic solution via single integrals
 * include cycle diagram(=>DHLL); (k0); (k2);(kn+n) 
 * =0, k-bond, =1, modified q-bond,etc, m>3,two diagrams are added up to evaluate
 * s = parameter * c^2 * I
 * ds/dc = parameter * ( 2*c*I + c^2 * dI/dc)
 * all are in simulation units
 * contains "linear step in sqrt(c)" && "linear step in c"
 * 
 * @author shu
 * date: 05/16/2012 
 */
public class B2Ionic {
	public static void main(String[] args) {
		double e = Electron.UNIT.toSim(1);
		double lambda = 4 * Math.PI/ Kelvin.UNIT.toSim(temperature) * e * e / relativePermitivity ; 
		System.out.println("2nd order diagrams");
		System.out.println("lambda at "+temperature+" K: "+lambda);
		System.out.println("sigmaHS:"+sigmaHS+"angstrom");
		double r_max = 600;// Defines range of separation distance, r = [0 r_max)
		int power = 16;
		int N = 1<<power;
		double del_r = r_max/N;
		double[] sqrt_c = new double[50];
		double[] c=new double[sqrt_c.length];
		for (int t=0; t<sqrt_c.length; t++) {
			
			if (increase_in_sqrt_c){
				//*********************** linear increase in sqrt(c)*********************
				System.out.println("linear increase in square root of concentration");
				sqrt_c[t] = (t+1.0) / 1000;
				c[t] = sqrt_c[t] * sqrt_c[t];// get c from sqrt_c
			} 
			else {
				//*********************** linear increase in concentration **************
				System.out.println("linear increase in concentration");
				c[t]=(t+1.0)/20000;
				sqrt_c[t] = Math.sqrt(c[t]);// get sqrt_c from c
			}

			double kappa = sqrt_c[t] * Math.sqrt(lambda);
			double d_kappa_dc = kappa/2/c[t];
			//System.out.println(sqrt_c[t]/Math.sqrt(2*0.000602214179)); // x-axis, square root of concentration in mol/L
			//System.out.println(c[t]);// concentration in sim unit
			//System.out.println(c[t]/0.000602214179+"mol/L, "+c[t]+"atom/angstrom^3"+"  kappa:"+ kappa);

			//********************* Cycle diagram with analytical solutions (Debye Huckel Limiting Law) ******************//
			double s_cycle =  kappa * kappa * kappa  / 12.0 / Math.PI;
			double ds_cycle_dkappa = kappa * kappa / 4.0 / Math.PI; 
			double ds_cycle = ds_cycle_dkappa * d_kappa_dc;
			double betaPressure_cycle =  s_cycle - c[t] * ds_cycle; 
			double phiMinus1_cycle = betaPressure_cycle/c[t];
			// double pressure =  - kappa * kappa * kappa / 24 / Math.PI;
		
			//System.out.println(s_cycle);
			//System.out.println(-ds_cycle);
			//System.out.println(betaPressure_cycle);
			//System.out.println(phiMinus1_cycle);

			
			//********************* k0 diagram using existing equations 2nd for HS ***************************************************//
			double s_k  =  -2.0 * Math.PI / 3 * sigmaHS * sigmaHS * sigmaHS * c[t] * c[t] ;
			double ds_k =  -4.0 * Math.PI / 3 * sigmaHS * sigmaHS * sigmaHS * c[t] ;
			double betaPressure_k =  s_k - c[t] * ds_k; 
			double phiMinus1_k = betaPressure_k/c[t];
			//System.out.println(phiMinus1_k);
			//System.out.println(-ds_k);
			//System.out.println(betaPressure_k);
			//System.out.println(s_k);

			
			//********************* k2 (one k-bond+two q-bonds) diagram using existing equations ***************************************************//
			double s_k2= kappa * kappa * kappa * ( Math.exp( -2 * kappa * sigmaHS) - 1 ) / 32 / Math.PI ; 
			double ds_k2 = ( 3 * kappa * kappa * ( Math.exp(-2 * kappa * sigmaHS)- 1 ) - 2 * sigmaHS * kappa * kappa * kappa * Math.exp( -2 * kappa * sigmaHS) )/ 32 / Math.PI * d_kappa_dc ;
			double betaPressure_k2 =  s_k2 - c[t] * ds_k2; 
			double phiMinus1_k2 = betaPressure_k2/c[t];
			//System.out.println(-ds_k2);
			//System.out.println(betaPressure_k2);
			//System.out.println(phiMinus1_k2);
			//System.out.println(s_k2);

			//*********************  q-bonds>2 diagram ***************************************************//
		    double[] fr  = getfr(N,del_r,kappa,c[t],bond,false);
			double[] dfr = getfr(N,del_r,kappa,c[t],bond,true);
			double parameter = Math.pow(lambda, bond) / 2 ;
			double sum = 0;
			double dsum = 0;
			
			// get integrals
			for (int i=0; i<fr.length; i++) {
				sum += fr[i]*del_r ;
				dsum += dfr[i]*del_r ;
			}
			double s = parameter * sum * c[t] * c[t] ;
			double ds=  parameter * ( 2 * c[t] * sum + c[t] * c[t] * dsum);
			double betaPressure= s - c[t] * ds;
			double phiMinus1 = betaPressure/c[t];
			//System.out.println(s);
			//System.out.println(-ds);
			//System.out.println(betaPressure);
			//System.out.println(phiMinus1);

			
		}
	}

	//Returns a discretized qm-like function except the domain [sigmaHS, infinity]
	public static double[] getfr(int N,double del_r, double kappa,double concentration,int bondNumber,boolean derivative) {
		double[] fr = new double[N]; 
		for (int n = 0; n<N; n++) {
			double r = n*del_r;
			if ( r > sigmaHS ) {
				fr[n] = 4 * Math.PI * r * r * Math.pow( ( Math.exp(-kappa * r)/ 4 / Math.PI / r ), bondNumber ) / SpecialFunctions.factorial(bondNumber); 
			}
			else {
				fr[n]=0;
			}
			if (derivative){
				fr[n]  *= -bondNumber * r * (kappa / 2 / concentration)  ; //dkappadc = kappa / 2 / c ;
			}
		}
		return fr;
	}
	private static double sigmaHS = 4.25 ;
	private static double temperature = 298;
	private static double relativePermitivity = 78.5 ; 
	private static int bond = 8; 
	private static boolean increase_in_sqrt_c = true;
}
