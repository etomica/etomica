/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.api.ISpecies;
import etomica.atom.IAtomList;
import etomica.space.Space;
import etomica.space.Vector;

import java.io.Serializable;

/**
 * Calculates the geometric center over a set of atoms. The position of the
 * atom or child atoms are accumulated and used to compute their
 * center (unweighted by mass). Calculated center is obtained via the getPosition
 * method.
 * This class is for normal alkane, TraPPE-Explicit Hydrogen only
 * Only carbons are employed to calculate the geometric center
 *
 * @author shu
 * March 2013
 */
public class MoleculePositionGeometricCenterAlkaneEH implements IMoleculePositionDefinition, Serializable {

    public MoleculePositionGeometricCenterAlkaneEH(Space space, ISpecies speciesAlkane) {
        center = space.makeVector();
        this.speciesAlkane =speciesAlkane; 
    }

    public Vector position(IMolecule molecule) {
        center.E(0.0);
        IAtomList children = molecule.getChildList();
        int nAtoms = children.getAtomCount();
        // get the species info in the virial main class, if it is alkaneEH species, then use (n-2)/3, use nAtoms for CO2/N2/etc.
        int numCarbons = molecule.getType()== speciesAlkane? (nAtoms-2)/3 : nAtoms;
        for (int i=0; i<numCarbons; i++) {// loop over all carbons ONLY
            center.PE(children.getAtom(i).getPosition());
        }
        center.TE(1.0 / numCarbons);
        return center;
    }

    private static final long serialVersionUID = 1L;
    private final Vector center;
    private final ISpecies speciesAlkane;
}
