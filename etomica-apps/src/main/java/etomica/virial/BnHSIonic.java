/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.virial.cluster.Standard;

/**
 * ionic solution, restrictive primitive model (RPM)
 * compute s, (phi[osmotic coefficient]-1) and ln(gamma) from HARD SPHERE interaction(from Standard class)
 * 2nd ~ 5th order integrals are included
 * p[] (<==> HSB[]) ==> integral = HSB[]* n! / (1-n)
 * s[n] = c^n * HSB[n] / (1-n)
 * ds[n]/dc = n * c ^ (n-1) * HSB[n] / (1-n) = s[n] * n/c  cuz integral is independent of concentration
 * in simulation units
 *  
 *  contains "linear step in sqrt(c)" && "linear step in c"
 *  
 * @author shu
 * date:April 27th 2012 
 */
public class BnHSIonic {
	public static void main(String[] args) {
		System.out.println("Ionic solution, RPM, contribution from hard sphere interactions,"+i+" order");
		System.out.println("ionic particle diameter: "+sigmaHS+"angstrom");
	
		double[] sqrt_c = new double[50];
		double[] c=new double[sqrt_c.length];
		
		for (int t=0; t<sqrt_c.length; t++) {
			
			//*********************** linear increase in sqrt(c)*********************
			//System.out.println("linear increase in square root of concentration");
			sqrt_c[t] = (t+1.0) / 1000;
			c[t] = sqrt_c[t] * sqrt_c[t];// get c from sqrt_c
				
			//*********************** linear increase in concentration **************
			//System.out.println("linear increase in concentration");
			//c[t]=(t+1.0)/20000;
			//sqrt_c[t] = Math.sqrt(c[t]);// get sqrt_c from c
			
			//*********************** the following is the same for both sqrt_c & c linear increase **************
			//System.out.println(c[t]);
			//System.out.println(c[t]/0.000602214179+"mol/L");
			//System.out.println(sqrt_c[t]/Math.sqrt(2*0.000602214179)); // x-axis, square root of c in mol/L
			
			double[] HSB = new double[6];//virial coefficient from standard class
			double[] s = new double[HSB.length];
			double[] ds = new double[HSB.length];
			double[] lnGamma = new double[HSB.length];
			double[] betaP = new double[HSB.length];
			double[] phiMinusOne = new double[HSB.length];
			HSB[2] = Standard.B2HS(sigmaHS);
			HSB[3] = Standard.B3HS(sigmaHS);
			HSB[4] = Standard.B4HS(sigmaHS);
			HSB[5] = Standard.B5HS(sigmaHS);
			if ( i > HSB.length-1){ // length is 6, so HSB5 is the highest order that can be handled
				throw new RuntimeException("I can handle up to 5th order only!");
	        } 
			s[i] = HSB[i] / (1-i) *  Math.pow(c[t], i) ;
			ds[i] = s[i] * i / c[t];
			lnGamma[i] = -ds[i];
			betaP[i] = s[i] - ds[i] * c[t];
			phiMinusOne[i] = betaP[i]/c[t];
			System.out.println(c[t]*sigmaHS*sigmaHS*sigmaHS+" "+c[t]+" 1/A^3, "+c[t]/0.000602214179+" mol/L,  lnGamma:  "+lnGamma[i]);
			//System.out.println(-ds[i]);
			//System.out.println(s[i]);
			//System.out.println(phiMinusOne[i]);
			//System.out.println("(phi-1)["+i+"]:"+phiMinusOne[i]);
		}			
	}
	private static double sigmaHS = 4.25;
	private static int i = 3;
}
