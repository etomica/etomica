/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.molecule.IMolecule;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Length;

/**
 * Meter for tabulation of the radius of gyration of a set of chain molecules. 
 * 
 * @author David Kofke
 */
public class MeterRadiusGyration extends DataSourceScalar {

    public MeterRadiusGyration(Space space) {
        super("Radius of Gyration", Length.DIMENSION);
        cm = space.makeVector();
        realPos = space.makeVector();
        dr = space.makeVector();
    }

    public double getDataAsScalar() {
        if (box == null)
            throw new IllegalStateException(
                    "must call setBox before using meter");
        Boundary boundary = box.getBoundary();
        int nLeafAtomsTot = 0;
        double r2Tot = 0.0;
        for (IMolecule molecule : box.getMoleculeList()) {
            // loop over molecules
            IAtomList childList = molecule.getChildList();
            if (childList.size() < 2) {
                // a monatomic molecule
                continue;
            }

            // find center of mass
            //do the first iterate explicitly, assume there is at least
            // one leaf atom
            IAtom firstAtom = childList.get(0);
            int nLeafAtoms = 1;
            realPos.E(firstAtom.getPosition());
            cm.E(realPos);
            Vector prevPosition = firstAtom.getPosition();
            for (int iChild = 1; iChild < childList.size(); iChild++) {
                IAtom a = childList.get(iChild);
                nLeafAtoms++;
                Vector position = a.getPosition();
                dr.Ev1Mv2(position, prevPosition);
                //molecule might be wrapped around the box.  calculate
                //the real difference in position
                boundary.nearestImage(dr);
                //realPos is now the effective position of a
                realPos.PE(dr);
                cm.PE(realPos);
                prevPosition = position;
            }
            cm.TE(1.0 / nLeafAtoms);
            // calculate Rg^2 for this chain
            double r2 = 0.0;
            prevPosition = firstAtom.getPosition();
            realPos.E(firstAtom.getPosition());
            for (int iChild = 0; iChild < childList.size(); iChild++) {
                IAtom a = childList.get(iChild);
                Vector position = a.getPosition();
                dr.Ev1Mv2(position, prevPosition);
                //molecule might be wrapped around the box.  calculate
                //the real difference in position
                boundary.nearestImage(dr);
                //realPos is now the effective position of a
                realPos.PE(dr);
                dr.Ev1Mv2(realPos, cm);// = realPos.M(cm);
                r2 += dr.squared();
                prevPosition = position;
            }
            r2Tot += r2;
            nLeafAtomsTot += nLeafAtoms;
        }
        return r2Tot / nLeafAtomsTot;
    }

    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }

    /**
     * @param box
     *            The box to set.
     */
    public void setBox(Box box) {
        this.box = box;
    }

    private Box box;
    private final Vector cm, realPos;
    private final Vector dr;

}
