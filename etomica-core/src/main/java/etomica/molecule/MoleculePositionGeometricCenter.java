/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.atom.IAtomList;
import etomica.space.Space;
import etomica.space.Vector;

import java.io.Serializable;

/**
 * Calculates the geometric center over a set of atoms. The position of the
 * atom or child atoms are accumulated and used to compute their
 * center (unweighted by mass). Calculated center is obtained via the getPosition
 * method.
 * 
 * @author David Kofke
 */
public class MoleculePositionGeometricCenter implements IMoleculePositionDefinition, Serializable {

    public MoleculePositionGeometricCenter(Space space) {
        center = space.makeVector();
    }

    public Vector position(IMolecule molecule) {
        center.E(0.0);
        IAtomList children = molecule.getChildList();
        int nAtoms = children.getAtomCount();
        for (int i=0; i<nAtoms; i++) {
            center.PE(children.getAtom(i).getPosition());
        }
        center.TE(1.0 / nAtoms);
        return center;
    }

    private static final long serialVersionUID = 1L;
    private final Vector center;
}
