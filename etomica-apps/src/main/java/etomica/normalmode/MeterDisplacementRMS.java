/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.api.IAtomList;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.data.DataSourceScalar;
import etomica.space.ISpace;
import etomica.units.Length;

public class MeterDisplacementRMS extends DataSourceScalar {

    protected final CoordinateDefinition coordinateDefinition;
    protected final IVectorMutable dr, boxSize0;
    
    public MeterDisplacementRMS(ISpace space, CoordinateDefinition coordinateDefinition) {
        super("displacement", Length.DIMENSION);
        this.coordinateDefinition = coordinateDefinition;
        dr = space.makeVector();
        boxSize0 = space.makeVector();
        boxSize0.E(coordinateDefinition.getBox().getBoundary().getBoxSize());
    }
    
    public double getDataAsScalar() {
        IAtomList leafList = coordinateDefinition.getBox().getLeafList();
        double sum = 0;
        IVector boxSize = coordinateDefinition.getBox().getBoundary().getBoxSize();
        for (int i=0; i<leafList.getAtomCount(); i++) {
            dr.E(coordinateDefinition.getLatticePosition(leafList.getAtom(i)));
            dr.DE(boxSize0);
            dr.TE(boxSize);
            dr.ME(leafList.getAtom(i).getPosition());
            sum += dr.squared();
//            if (i==0) System.out.println("hi there "+Math.sqrt(dr.squared()));
        }
        return Math.sqrt(sum/leafList.getAtomCount()/3);
    }
}
