/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential.amoeba;

import etomica.potential.IPotential2;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;

/**
 * Spherically symmetric potential of the form u(r) = Kb*(r-r0)^2*(1 - 2.55 * (r - r0) + 7/12 * 2.55 * (r-r0)^2)
 */
public class P2AmoebaStretch implements IPotential2 {

    private double Kb;// Spring constant gives a measure of the strength of harmonic interaction
	private double r0;

    /**
     *
     * @param Kb energy constant
     * @param r0  Separation at which potential is at its minimum.  Default is
     * zero.
     */
    public P2AmoebaStretch(double Kb, double r0) {
        setK(Kb);
        setR0(r0);
    }

    public void u012add(double r2, double[] u012) {
        double r = Math.sqrt(r2);
        double dx = r - r0;
        u012[0] = Kb*dx*dx*(1 - 2.55 * dx + 7.0/12.0*2.55*2.55 * dx*dx);
        u012[1] = Kb*r*dx*(2 - 3*2.55*dx + 4*7.0/12.0*2.55*2.55 * dx*dx);
//        Unit kcalpmole = new UnitRatio(new PrefixedUnit(Prefix.KILO, Calorie.UNIT), Mole.UNIT);
//        System.out.println("bond "+Math.sqrt(r2)+" "+r0+" "+dx+" "+kcalpmole.fromSim(Kb)+" "+kcalpmole.fromSim(u012[0]));
    }

    public double u(double r2) {
    	double dx = Math.sqrt(r2) - r0;
    	return Kb*dx*dx*(1 - 2.55 * dx + 7.0/12.0*2.55 * dx*dx);
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
    	double r = Math.sqrt(r2);
        double dx = r - r0;
    	return Kb*r*dx*(2 - 3*2.55*dx + 4*7.0/12.0*2.55 * dx*dx);
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
        throw new RuntimeException("Nope");
    }

    /**
     * Accessor method for harmonic energy parameter
     */
    public double getK() {return Kb;}
    /**
     * Accessor method for harmonic energy parameter
     */
    public void setK(double K) {
        Kb = K;
    }

    /**
     * Not implemented correctly.  
     * Should be energy/length^2.
     */
    public Dimension getSpringConstantDimension() {
        return new CompoundDimension(new Dimension[]{Energy.DIMENSION,Length.DIMENSION},new double[]{1,-2});
    }
    
	/**
	 * Separation at which potential is at its minimum.
	 * @return double
	 */
	public double getR0() {
		return r0;
	}

	/**
	 * Sets the the separation at which potential is at its minimum.
	 * @param r0 The r0 to set
	 */
	public void setR0(double r0) {
		this.r0 = r0;
	}

    public Dimension getR0Dimension() {
        return Length.DIMENSION;
    }
    
}
