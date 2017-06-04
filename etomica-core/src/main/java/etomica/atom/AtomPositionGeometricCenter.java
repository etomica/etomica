/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import java.io.Serializable;

import etomica.api.IMolecule;
import etomica.space.Vector;
import etomica.space.Space;

/**
 * Calculates the geometric center over a set of atoms. The position of the
 * atom or child atoms are accumulated and used to compute their
 * center (unweighted by mass). Calculated center is obtained via the getPosition
 * method.
 * 
 * @author David Kofke
 */
public class AtomPositionGeometricCenter implements IAtomPositionDefinition, Serializable {

    public AtomPositionGeometricCenter(Space space) {
        center = space.makeVector();
    }

    public Vector position(IMolecule atom) {
        center.E(0.0);
        IAtomList children = atom.getChildList();
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
