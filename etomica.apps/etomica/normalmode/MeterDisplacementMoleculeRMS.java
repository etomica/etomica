/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.api.IMoleculeList;
import etomica.api.IVectorMutable;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.IAtomPositionDefinition;
import etomica.data.DataSourceScalar;
import etomica.space.ISpace;
import etomica.units.Length;

public class MeterDisplacementMoleculeRMS extends DataSourceScalar {

    protected final CoordinateDefinitionMolecule coordinateDefinition;
    protected final IVectorMutable dr;
    protected IAtomPositionDefinition position;
    
    public MeterDisplacementMoleculeRMS(ISpace space, CoordinateDefinitionMolecule coordinateDefinition) {
        super("displacement", Length.DIMENSION);
        this.coordinateDefinition = coordinateDefinition;
        dr = space.makeVector();
        position = new AtomPositionGeometricCenter(space);
    }
    
    public double getDataAsScalar() {
        IMoleculeList molecules = coordinateDefinition.getBox().getMoleculeList();
        double sum = 0;
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            dr.E(position.position(molecules.getMolecule(i)));
            dr.ME(coordinateDefinition.getLatticePosition(molecules.getMolecule(i)));
            sum += dr.squared();
        }
        return Math.sqrt(sum/molecules.getMoleculeCount()/3);
    }
}
