/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.materialfracture;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.units.Pressure2D;

/**
 * Meter that calculates the stress within the gage cell.
 *
 * @author Andrew Schultz
 */
public class MeterLoad extends DataSourceScalar {
    
    public MeterLoad(PotentialCalculationForceStress pc) {
        super("Load", Pressure2D.DIMENSION);
        this.pc = pc;
    }

    public double getDataAsScalar(){
        return pc.getLoad();
    }

    protected final PotentialCalculationForceStress pc;
    protected Box box;
}
