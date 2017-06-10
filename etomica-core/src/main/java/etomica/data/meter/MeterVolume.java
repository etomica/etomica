/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.units.dimensions.Volume;

/**
 * Meter for measurement of the box volume.
 */
public class MeterVolume extends DataSourceScalar {
    
    public MeterVolume() {
        super("Volume",Volume.DIMENSION);
    }

    public double getDataAsScalar() {
        return box.getBoundary().volume();
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

    private static final long serialVersionUID = 1L;
    private Box box;
}
