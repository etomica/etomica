/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.rheology;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;

/**
 * Meter that measures that end to end distance of the molecules
 *
 * @author Andrew Schultz
 */
public class MeterEndToEnd extends DataSourceScalar {

    public MeterEndToEnd(Space space) {
        super("end-to-end distance^2", new CompoundDimension(new Dimension[]{Length.DIMENSION}, new double[]{2}));
        dr = space.makeVector();
    }

    public void setBox(Box newBox) {
        box = newBox;
    }
    
    public double getDataAsScalar() {
        IMoleculeList molecules = box.getMoleculeList();
        double ee_tot = 0;
        for (int i = 0; i<molecules.size(); i++) {
            IAtomList atoms = molecules.get(i).getChildList();
            dr.E(atoms.get(atoms.size()-1).getPosition());
            dr.ME(atoms.get(0).getPosition());
            box.getBoundary().nearestImage(dr);
            ee_tot += dr.squared();
        }
        return ee_tot/molecules.size();
    }

    protected Box box;
    protected final Vector dr;
}
