/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.units.dimensions.Energy;

/**
 * Meter that returns the potential energy using an IAPIPotential.
 * 
 * @author Andrew Schultz
 */
public class MeterAPIPotentialEnergy extends DataSourceScalar {
    
    public MeterAPIPotentialEnergy(IAPIPotential potential) {
        super("energy", Energy.DIMENSION);
        this.potential = potential;
    }

    public double getDataAsScalar() {
        return potential.calculateEnergy(box);
    }
    
    public void setBox(Box newBox) {
        box = newBox;
    }
    
    public Box getBox() {
        return box;
    }

    private static final long serialVersionUID = 1L;
    protected final IAPIPotential potential;
    protected Box box;
}
