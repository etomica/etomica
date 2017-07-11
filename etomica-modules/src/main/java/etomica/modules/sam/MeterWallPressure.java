/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.sam;

import etomica.box.Box;
import etomica.space.Vector;
import etomica.data.DataSourceScalar;
import etomica.units.dimensions.Pressure;

public class MeterWallPressure extends DataSourceScalar {


    public MeterWallPressure(PotentialCalculationForceSumWall pc) {
        super("Wall Pressure", Pressure.DIMENSION);
        this.pc = pc;
    }

    public void setBox(Box newBox) {
        box = newBox;
    }
    
    public double getDataAsScalar() {
        double f = pc.getWallForce();
        Vector dimensions = box.getBoundary().getBoxSize();
        double A = 1;
        for (int i=0; i<dimensions.getD(); i++) {
            if (i == pc.getWallPotential().getWallDim()) {
                continue;
            }
            A *= dimensions.getX(i);
        }
        return f / A;
    }

    private static final long serialVersionUID = 1L;
    protected final PotentialCalculationForceSumWall pc;
    protected Box box;
}
