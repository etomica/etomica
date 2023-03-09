/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.amoeba;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotential2;
import etomica.potential.IPotentialMolecular;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesManager;
import etomica.util.collections.IntArrayList;

public class PotentialMoleculePairAmoebaVDW extends PotentialMoleculeAmoebaVDW implements IPotentialMolecular {

    public PotentialMoleculePairAmoebaVDW(Space space, SpeciesManager sm, IntArrayList[][] bonding) {
        super(space, sm, bonding);
    }

    public PotentialMoleculePairAmoebaVDW(Space space, IPotential2[][] atomPotentials, IntArrayList[][] bonding) {
        super(space, atomPotentials, bonding);
    }

    @Override
    public double energy(IMoleculeList molecules) {
        return energy(molecules.get(0), molecules.get(1));
    }


    public double energy(IMolecule molecule1, IMolecule molecule2) {
        IAtomList atoms1 = molecule1.getChildList();
        IAtomList atoms2 = molecule2.getChildList();
        double u = 0;
        Vector[] r2 = new Vector[atoms2.size()];
        for (IAtom a2 : atoms2) {
            r2[a2.getIndex()] = getReducedPosition(a2);
        }
        for (IAtom a1 : atoms1) {
            Vector r1 = getReducedPosition(a1);
            IPotential2[] p1 = atomPotentials[a1.getType().getIndex()];
            for (IAtom a2 : atoms2) {
                IPotential2 p2 = p1[a2.getType().getIndex()];
                if (p2 == null) continue;
                Vector dr = space.makeVector();
                dr.Ev1Mv2(r2[a2.getIndex()], r1);
                double uu = p2.u(dr.squared());
                u += uu;
            }
        }
        return u;
    }

}
