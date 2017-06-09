/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.api.IPotentialAtomic;
import etomica.api.IPotentialMolecular;
import etomica.atom.IAtomList;
import etomica.molecule.IMoleculeList;

import java.math.BigDecimal;
import java.math.MathContext;

/**
 * Mimics original PotentialEnergySum calculation, but uses a BigDecimal
 * representation for the sum.  This can eliminate numerical roundoff error
 * caused by summing over molecule pairs.
 * 
 * @author Andrew Schultz
 */
public class PotentialCalculationEnergySumBigDecimal extends PotentialCalculationEnergySum {

    public PotentialCalculationEnergySumBigDecimal(int precision) {
        sum = BigDecimal.ZERO;
        this.precision = precision;
        mc = new MathContext(precision);
    }
    
    /**
	 * Adds to the energy sum the energy values obtained from application of the given potential to the
	 * atoms produced by the given iterator.  Iterator is reset by method before beginning calculation.
	 */
	public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
	    BigDecimal u = new BigDecimal(potential.energy(atoms), mc);
	    sum = sum.add(u); 
	}
	
    /**
     * Adds to the energy sum the energy values obtained from application of the given potential to the
     * atoms produced by the given iterator.  Iterator is reset by method before beginning calculation.
     */
    public void doCalculation(IMoleculeList atoms, IPotentialMolecular potential) {
        BigDecimal u = new BigDecimal(potential.energy(atoms), mc);
        sum = sum.add(u); 
    }
    
	/**
	 * Sets the energy sum to zero, typically to begin a new energy-sum calculation.
	 */
	public void zeroSum() {
		sum = BigDecimal.ZERO;
	}

	/**
	 * Returns the current value of the energy sum.
	 */
	public double getSum() {
        return sum.doubleValue();
    }
	
	public int getPrecision() {
	    return precision;
	}
	
	public void setPrecision(int newPrecision) {
	    precision = newPrecision;
	    mc = new MathContext(newPrecision);
	}
	
    private static final long serialVersionUID = 1L;
	protected BigDecimal sum;
	protected int precision;
	protected MathContext mc;
}
