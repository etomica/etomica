/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;

/**
 * Calculates the center of mass (COM) for a molecule in a box.  The class
 * handles molecules that span the box by applying nearest image to each atom's
 * vector displacement from the previously-computed center.  This should work
 * so long as the molecule doesn't wrap back on itself and will always work if
 * the molecules extent is less than L/2.
 */
public class CenterOfMass {

    public static Vector position(Box box, IMolecule molecule) {
        double massSum = 0;
        Vector center = box.getSpace().makeVector();
        Vector dr = Vector.d(center.getD());
        IAtomList children = molecule.getChildList();
        int nAtoms = children.size();
        for (int i = 0; i < nAtoms; i++) {
            IAtom a = children.get(i);
            double mass = a.getType().getMass();
            if (massSum == 0) {
                center.PEa1Tv1(mass, a.getPosition());
            } else {
                // sum = sum + mass*((sum/n)+pbc(r - sum/n))
                dr.E(a.getPosition());
                dr.PEa1Tv1(-1 / massSum, center);
                box.getBoundary().nearestImage(dr);
                dr.PEa1Tv1(1 / massSum, center);
                center.PEa1Tv1(mass, dr);
            }
            massSum += mass;
        }
        center.TE(1.0 / massSum);
        center.PE(box.getBoundary().centralImage(center));
        return center;
    }
}
