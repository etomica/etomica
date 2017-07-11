/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.atom.iterator.AtomIteratorBoxDependent;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.DataSourceScalar;
import etomica.math.geometry.Polytope;
import etomica.space.Boundary;
import etomica.units.dimensions.DimensionRatio;
import etomica.units.dimensions.Quantity;
import etomica.units.dimensions.Volume;

/**
 * Meter for measurement of density within a specified subvolume
 */
 
public abstract class MeterLocalDensity extends DataSourceScalar {

    public MeterLocalDensity(Box _box) {
        super("Local Density",new DimensionRatio(Quantity.DIMENSION, Volume.DIMENSION));
        box = _box;
        iterator = new AtomIteratorLeafAtoms(box);
        if (box.getBoundary() instanceof Boundary) {
            setShape(((Boundary)box.getBoundary()).getShape());
        }
    }

    public void setShape(Polytope shape) {
        this.shape = shape;
    }

    public Polytope getShape() {
        return shape;
    }

    /**
     * @return the current value of the local density or local mole fraction
     */
    public double getDataAsScalar() {
        if (shape == null) throw new IllegalStateException("must call setShape before using meter (or use a Boundary)");
        //compute local molar density
        int nSum = 0;
        iterator.reset();
        for (IAtom atom = iterator.nextAtom(); atom != null;
             atom = iterator.nextAtom()) {
            if(shape.contains(atom.getPosition())) nSum++;
        }
        return nSum/shape.getVolume();
    }

    /**
     * @return Returns the iterator.
     */
    public AtomIteratorBoxDependent getIterator() {
        return iterator;
    }
    /**
     * @param iterator The iterator to set.
     */
    public void setIterator(AtomIteratorBoxDependent iterator) {
        this.iterator = iterator;
    }

    private static final long serialVersionUID = 1L;
    private Box box;
    /**
     * Class variable used to specify that all species are included in number-density calculation
     */
    private AtomIteratorBoxDependent iterator;
    private Polytope shape;
}
