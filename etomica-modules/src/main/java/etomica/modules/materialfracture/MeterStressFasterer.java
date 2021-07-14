/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.materialfracture;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.space.Vector;
import etomica.units.dimensions.Pressure2D;

/**
 * Meter that calculates the stress within the gage cell.
 *
 * @author Andrew Schultz
 */
public class MeterStressFasterer extends DataSourceScalar {

    public MeterStressFasterer(MaterialFractureFasterer sim) {
        super("Stress", Pressure2D.DIMENSION);
        this.sim = sim;
    }

    public void setBox(Box newBox) {
        box = newBox;
    }

    public Box getBox() {
        return box;
    }

    public double getDataAsScalar() {
        double area = 1;
        Vector dim = box.getBoundary().getBoxSize();
        for (int i = 1; i < dim.getD(); i++) {
            area *= dim.getX(i);
        }

        return sim.wallForce / area / 2.0;
    }

    protected final MaterialFractureFasterer sim;
    protected Box box;
}
