/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IMoleculeList;
import etomica.atom.IMoleculePositionDefinition;
import etomica.space.Vector;
import etomica.atom.MoleculePositionGeometricCenter;
import etomica.data.DataSourceScalar;
import etomica.space.Space;
import etomica.units.Length;

public class MeterDisplacementMoleculeRMS extends DataSourceScalar {

    protected final CoordinateDefinitionMolecule coordinateDefinition;
    protected final Vector dr;
    protected IMoleculePositionDefinition position;
    
    public MeterDisplacementMoleculeRMS(Space space, CoordinateDefinitionMolecule coordinateDefinition) {
        super("displacement", Length.DIMENSION);
        this.coordinateDefinition = coordinateDefinition;
        dr = space.makeVector();
        position = new MoleculePositionGeometricCenter(space);
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
