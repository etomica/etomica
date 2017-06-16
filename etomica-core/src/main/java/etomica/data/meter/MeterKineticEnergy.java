/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.atom.IAtomKinetic;
import etomica.box.Box;
import etomica.atom.iterator.AtomIteratorBoxDependent;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.DataSourceScalar;
import etomica.units.Energy;

/**
 * Meter for the total kinetic energy in a box
 * Computes total KE by summing values of KE returned by every atom in the box.
 * A different box-dependent atom integrator may be set to permit calculation
 * over a particular set of atoms in the box.
 */
 
public class MeterKineticEnergy extends DataSourceScalar {
    private static final long serialVersionUID = 1L;
    private AtomIteratorBoxDependent iterator;
    
    public MeterKineticEnergy() {
        super("Kinetic Energy",Energy.DIMENSION);
        setIterator(new AtomIteratorLeafAtoms());
    }
    
    /**
     * Returns the iterator that defines the atoms summed for their
     * kinetic energy.
     */
	public AtomIteratorBoxDependent getIterator() {
		return iterator;
	}
	
	/**
	 * Sets the iterator that defines the atoms which are summed for
	 * their total kinetic energy.  Default is a leaf-atom iterator,
	 * giving all leaf atoms in the box.
	 * @param iterator
	 */
	public void setIterator(AtomIteratorBoxDependent iterator) {
		this.iterator = iterator;
	}
	
	/**
	 * Returns the total kinetic energy summed over all atoms produced by
	 * the iterator when applied to the given box.  Does not include contributions
     * from atoms having infinite mass (it assumes they are stationary).
	 */
    public double getDataAsScalar() {
        if (box == null) throw new IllegalStateException("must call setBox before using meter");
        double ke = 0.0;
        iterator.setBox(box);
        iterator.reset();
        for (IAtomKinetic atom = (IAtomKinetic)iterator.nextAtom(); atom != null;
             atom = (IAtomKinetic)iterator.nextAtom()) {
            double mass = atom.getType().getMass();
            if(mass == Double.POSITIVE_INFINITY) continue;
            ke += 0.5*mass*(atom.getVelocity().squared());
        }
        return ke;
    }
    
    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }
    /**
     * @param box The box to set.
     */
    public void setBox(Box box) {
        this.box = box;
    }

    private Box box;
 }
