/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.molecule.IMoleculeList;
import etomica.space.Vector;

/**
 * Evaluates the energy summed over all iterated atoms. Each call to doCalculate
 * accumulates the energy of the given potential applied to the atoms produced
 * by the given iterator. Each such call accumulates a single sum, which can
 * be accessed via the getSum() method.  Sum is re-zeroed only upon a call to the
 * zeroSum() method. 
 *
 * @author David Kofke
 */
public class PotentialCalculationEnergySum implements PotentialCalculation, PotentialCalculationMolecular, java.io.Serializable {

    public static boolean debug = false;
    
    /**
	 * Adds to the energy sum the energy values obtained from application of the given potential to the
	 * atoms.
	 */
	public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
	    sum += potential.energy(atoms);
	    if (debug && Double.isInfinite(sum) || Double.isNaN(sum)) {
	        System.err.println("unhappy energy "+sum+" for "+atoms+" "+atoms.getAtom(0).hashCode()+" "+atoms.getAtom(1).hashCode()+" "+" "+atoms.getAtom(0).getPosition()+" "+atoms.getAtom(1).getPosition());
	        Vector v1 = atoms.getAtom(0).getPosition();
			Vector v2 = atoms.getAtom(1).getPosition();
			double distance  = Math.sqrt(v1.Mv1Squared(v2));
			System.err.println("distance "+distance);
	        potential.energy(atoms);
	        debug = false;
	    }
	}
	
    /**
     * Adds to the energy sum the energy values obtained from application of the given potential to the
     * molecules.
     */
    public void doCalculation(IMoleculeList molecules, IPotentialMolecular potential) {
        sum += potential.energy(molecules);
        if (debug && Double.isInfinite(sum) || Double.isNaN(sum)) {
            System.err.println("unhappy energy "+sum+" for "+molecules);
            potential.energy(molecules);
            debug = false;
        }
    }
    
	/**
	 * Sets the energy sum to zero, typically to begin a new energy-sum calculation.
	 */
	public void zeroSum() {
		sum = 0.0;
	}

	/**
	 * Returns the current value of the energy sum.
	 */
	public double getSum() {
        return sum;
    }
	
    private static final long serialVersionUID = 1L;
	protected double sum = 0.0;
}
