/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.math.geometry.Polytope;
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
        setShape(box.getBoundary().getShape());
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
        for (IAtom atom : box.getLeafList()) {
            if(shape.contains(atom.getPosition())) nSum++;
        }
        return nSum/shape.getVolume();
    }

    private Box box;
    /**
     * Class variable used to specify that all species are included in number-density calculation
     */
    private Polytope shape;
}
