/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Boundary;
import etomica.space.Vector;

/**
 * Calculates the center of mass (COM) over a set of atoms accounting for
 * periodic boundary conditions.  Molecules that wrap around the box will not
 * have a COM in the middle of the box and the returned COM will always be in
 * the box.
 */
public class MoleculePositionCOMPBC implements IMoleculePositionDefinition {

    protected final Boundary boundary;

    public MoleculePositionCOMPBC(Box box) {
        boundary = box.getBoundary();
    }

    public Vector position(IMolecule molecule) {
        return com(boundary, molecule);
    }

    public static Vector com(Boundary boundary, IMolecule molecule) {
        double massSum = 0;
        Vector center = null, dr = null, pos0 = null;
        IAtomList children = molecule.getChildList();
        int nAtoms = children.size();
        for (int i=0; i<nAtoms; i++) {
            IAtom a = children.get(i);
            double mass = a.getType().getMass();
            Vector p = a.getPosition();
            if (i==0) {
                pos0 = p;
                center = Vector.d(p.getD());
                dr = Vector.d(p.getD());
            }
            else {
                dr.Ev1Mv2(children.get(i).getPosition(), pos0);
                boundary.nearestImage(dr);
                center.PEa1Tv1(mass, dr);
            }
            massSum += mass;
        }
        center.TE(1.0 / massSum);
        center.PE(pos0);
        center.PE(boundary.centralImage(center));
        return center;
    }
}
