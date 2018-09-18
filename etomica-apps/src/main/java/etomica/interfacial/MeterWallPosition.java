/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.interfacial;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.species.ISpecies;
import etomica.units.dimensions.Length;

public class MeterWallPosition extends DataSourceScalar {

    protected final Box box;
    protected final ISpecies wallSpecies;
    protected final double zShift;
    
    public MeterWallPosition(Box box, ISpecies wallSpecies, double zShift) {
        super("Wall position", Length.DIMENSION);
        this.box = box;
        this.wallSpecies = wallSpecies;
        this.zShift = zShift;
    }

    public double getDataAsScalar() {
        return box.getMoleculeList(wallSpecies).get(0).getChildList().get(0).getPosition().getX(2) - zShift;
    }

}
