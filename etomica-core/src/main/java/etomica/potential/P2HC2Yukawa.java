/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Boundary;

/**
 * Hard-core plus two Yukawa fluid (HC2Yukawa): A Lennard-Jones like potential.
 * 
 * ==================================================================================================================================
 * 2001. Pini, Stell, and Wilding. "Liquid-gas box behaviour of an argon-like fluid modelled by the hard-core two-Yukawa potential"
 * ==================================================================================================================================
 * 	
 * 			| infinity																		r <= sigma
 * U(r) =	| 
 * 			|(A1/r) * exp[-z1 * (r - sigma)] - (A2 * epsilon/r) * exp[-z2 * (r - sigma)]	r > sigma
 *
 * 	LJ behavior parameters:	A1 = 1.6438 * sigma
 * 							z1 = 14.7 / sigma
 * 							A2 = 2.03 * sigma
 * 							z2 = 2.69 / sigma
 *
 * @author msellers
 */

public final class P2HC2Yukawa implements IPotential2 {

	public static IPotential2 makeTruncated(double sigma, double epsilon, TruncationFactory tf) {
		return tf.make(new P2HC2Yukawa(sigma, epsilon));
	}

	public double d2u(double r2) {
		// TODO Auto-generated method stub
		throw new RuntimeException();
	}

	public double du(double r2) {
		// TODO Auto-generated method stub
		throw new RuntimeException();
	}

	public P2HC2Yukawa(double sigma, double epsilon){
		setSigma(sigma);
		setEpsilon(epsilon);
		setParameters(sigma);
	}

	/**
	 * Energy method.  u(double r) form.
	 */
	public double u(double r2){
		double r = Math.sqrt(r2);
		//hard core repulsion
		if (r <= sigma){
			return Double.POSITIVE_INFINITY;
		}
		
		//two-Yukawa tail attraction
		double expRsigma = Math.exp(r - sigma);
		double uterm1 = (A1 / r) * expZ1 * expRsigma;
		double uterm2 = (A2 * epsilon / r) * expZ2 * expRsigma;
		
		return (uterm1 - uterm2);
	}

	/**
	 * Accessor methods for size and energy parameters.
	 */
	public double getSigma() {return sigma;}
	
	public double getEpsilon() {return epsilon;}
	
	/**
	 * Mutator methods for size and energy parameters.
	 */
	public final void setEpsilon(double eps) {epsilon = eps;}
	
	public final void setSigma(double s) {sigma = s;}
	
	/**
	 * Parameter calculation method.
	 */
	public final void setParameters(double s){
		A1 = 1.6438 * s;
		A2 = 2.03 * s;
		z1 = 14.7 / s;
		z2 = 2.69 / s;
		expZ1 = Math.exp(-z1);
		expZ2 = Math.exp(-z2);
	}
	
	private double sigma;
	private double epsilon;
	private double A1;
	private double A2;
	private double z1;
	private double z2;
	private double expZ1;
	private double expZ2; 
	private Boundary nearestImageTransformer;
}	
