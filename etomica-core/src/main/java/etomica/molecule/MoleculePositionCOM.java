/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Space;
import etomica.space.Vector;

import java.io.Serializable;

/**
 * Calculates the center of mass (COM) over a set of atoms.
 *
 * @author David Kofke
 */

public class MoleculePositionCOM implements IMoleculePositionDefinition, Serializable {

    public MoleculePositionCOM(Space space) {
        center = space.makeVector();
    }
    
    public Vector position(IMolecule molecule) {
        double massSum = 0;
        center.E(0.0);
        IAtomList children = molecule.getChildList();
        int nAtoms = children.size();
        for (int i=0; i<nAtoms; i++) {
            IAtom a = children.get(i);
            double mass = a.getType().getMass();
            center.PEa1Tv1(mass, a.getPosition());
            massSum += mass;
        }
        center.TE(1.0 / massSum);
        return center;
    }

    private static final long serialVersionUID = 1L;
    private final Vector center;
}
