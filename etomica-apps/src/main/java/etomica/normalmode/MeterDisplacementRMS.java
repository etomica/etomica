/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.data.DataSourceScalar;
import etomica.space.Space;
import etomica.units.dimensions.Length;

public class MeterDisplacementRMS extends DataSourceScalar {

    protected final CoordinateDefinition coordinateDefinition;
    protected final Vector dr, boxSize0;
    
    public MeterDisplacementRMS(Space space, CoordinateDefinition coordinateDefinition) {
        super("displacement", Length.DIMENSION);
        this.coordinateDefinition = coordinateDefinition;
        dr = space.makeVector();
        boxSize0 = space.makeVector();
        boxSize0.E(coordinateDefinition.getBox().getBoundary().getBoxSize());
    }
    
    public double getDataAsScalar() {
        IAtomList leafList = coordinateDefinition.getBox().getLeafList();
        double sum = 0;
        Vector boxSize = coordinateDefinition.getBox().getBoundary().getBoxSize();
        for (int i = 0; i<leafList.size(); i++) {
            dr.E(coordinateDefinition.getLatticePosition(leafList.get(i)));
            dr.DE(boxSize0);
            dr.TE(boxSize);
            dr.ME(leafList.get(i).getPosition());
            sum += dr.squared();
//            if (i==0) System.out.println("hi there "+Math.sqrt(dr.squared()));
        }
        return Math.sqrt(sum/leafList.size()/3);
    }
}
