/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import java.io.Serializable;

import etomica.space.Vector;
import etomica.space.Space;

/**
 * Calculates the center of mass (COM) over a set of atoms. The mass and
 * position of all atoms passed to the actionPerformed are accumulated and used
 * to compute their center of mass. Calculated COM is obtained via the getData
 * method, which returns an array with the elements of the calculated COM
 * vector, or via the getCOM method, which returns a vector with the calculated
 * COM. Calculation is zeroed via the reset method.
 * <p>
 * A typical use of this class would have it passed to the allAtoms method of an
 * iterator, or wrapped in an AtomGroupAction to calculate the COM of the atoms
 * in an atom group.
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
        int nAtoms = children.getAtomCount();
        for (int i=0; i<nAtoms; i++) {
            IAtom a = children.getAtom(i);
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
