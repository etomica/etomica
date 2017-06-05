/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.api.*;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.units.Angle;

/**
 * Meter that measures the average tilt angle (not the angle of average tilt!)
 *
 * @author Andrew Schultz
 */
public class MeterPhiDeviation extends DataSourceScalar {

    public MeterPhiDeviation(Space space) {
        super("phi deviation", Angle.DIMENSION);
        dr = space.makeVector();
    }
    
    public void setBox(Box newBox) {
        box = newBox;
    }

    public double getDataAsScalar() {
        IMoleculeList molecules = box.getMoleculeList();
        int nMolecules = molecules.getMoleculeCount();
        double sum = 0;
        for (int i=0; i<nMolecules; i++) {
            IMolecule molecule = molecules.getMolecule(i);
            IAtomList atomList = molecule.getChildList();
            int leafCount = atomList.getAtomCount();
            dr.E(atomList.getAtom(leafCount-1).getPosition());
            dr.ME(atomList.getAtom(0).getPosition());
            dr.normalize();
            double phi = Math.atan2(dr.getX(1), dr.getX(0));
            double sintheta = Math.sqrt(dr.getX(0)*dr.getX(0) + dr.getX(1)*dr.getX(1));
            double u = phi*sintheta;
            sum += u*u;
        }
        return Math.sqrt(sum/nMolecules);
    }

    private static final long serialVersionUID = 1L;
    protected Box box;
    protected final Vector dr;
}
