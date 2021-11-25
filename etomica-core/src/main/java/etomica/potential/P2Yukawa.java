/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Space;

/**
 * Yukawa interatomic potential.
 *
 * U(r) = Vo * exp(-K * r) / r
 *
 * where Vo is the potential's energy parameter and K is the characteristic distance (1/K is a measure for the screening length).
 *
 * @author msellers
 *
 */

public final class P2Yukawa implements Potential2Soft {

	public static Potential2Soft makeTruncated(double kappa, double vzero, TruncationFactory tf) {
		return tf.make(new P2Yukawa(kappa, vzero));
	}

	public P2Yukawa() {
		this(1.0, 1.0);
	}

	public P2Yukawa(double kappa, double vzero) {
		setKappa(kappa);
		setVZero(vzero);
	}

	/**
	 * Energy method.  u(double r) form.
	 */
	public double u(double r2){
		double r = Math.sqrt(r2);
		return vzero * Math.exp(-kappa * r) / r;
	}
	
	/**
	 * r * du/dr method.
	 */
	public double du(double r2){
		double r = Math.sqrt(r2);
		return (-vzero * Math.exp(-kappa * r)) * (kappa + (1 / r));
	}
	
	/**
	 * r^2 * d^2u/dr^2 method.
	 */
	public double d2u(double r2){
		double r = Math.sqrt(r2);
		return (vzero * Math.exp(-kappa * r) * r) * (kappa * (kappa + (2 / r)) + (2 / r2));
	}
	
	
    /**
     * Integral from rC to infinity.
     */
    public double integral(Space space, double rC) {

        return 0;
    }

	/**
	 * Accessor methods for size and energy parameters.
	 */
	public double getKappa() {return kappa;}
	
	public double getVZero() {return vzero;}
	
	/**
	 * Mutator methods for size and energy parameters.
	 */
	public final void setKappa(double k) {kappa = k;}
	
	public final void setVZero(double v) {vzero = v;}

	
    private static final long serialVersionUID = 1L;
	private double kappa;
	private double vzero;
}

