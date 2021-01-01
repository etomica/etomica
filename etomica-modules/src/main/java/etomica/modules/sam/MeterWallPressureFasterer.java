/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.sam;

import etomica.data.DataSourceScalar;
import etomica.space.Vector;
import etomica.units.dimensions.Pressure;

public class MeterWallPressureFasterer extends DataSourceScalar {

    protected final SamFasterer sam;

    public MeterWallPressureFasterer(SamFasterer sam) {
        super("Wall Pressure", Pressure.DIMENSION);
        this.sam = sam;
    }

    public double getDataAsScalar() {
        double f = sam.wallForce;
        Vector dimensions = sam.box.getBoundary().getBoxSize();
        double A = 1;
        for (int i = 0; i < dimensions.getD(); i++) {
            if (i == 1) {
                continue;
            }
            A *= dimensions.getX(i);
        }
        return f / A;
    }
}
