/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.amoeba;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.BondingInfo;
import etomica.potential.IPotential2;
import etomica.potential.IPotentialMolecular;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesManager;
import etomica.util.collections.IntArrayList;

public class PotentialMoleculeIntraAmoebaVDW extends PotentialMoleculeAmoebaVDW implements IPotentialMolecular {

    protected final BondingInfo bondingInfo;

    public PotentialMoleculeIntraAmoebaVDW(Space space, SpeciesManager sm, IntArrayList[][] bonding, BondingInfo bondingInfo) {
        super(space, sm, bonding);
        this.bondingInfo = bondingInfo;
    }

    public PotentialMoleculeIntraAmoebaVDW(Space space, IPotential2[][] atomPotentials, IntArrayList[][] bonding, BondingInfo bondingInfo) {
        super(space, atomPotentials, bonding);
        this.bondingInfo = bondingInfo;
    }

    @Override
    public double energy(IMoleculeList molecules) {
        return energy(molecules.get(0));
    }


    public double energy(IMolecule molecule) {
        IAtomList atoms = molecule.getChildList();
        double u = 0;
        Vector[] r = new Vector[atoms.size()];
        for (IAtom a : atoms) {
            r[a.getIndex()] = getReducedPosition(a);
        }
        for (int i=0; i<atoms.size(); i++) {
            IAtom a1 = atoms.get(i);
            Vector r1 = r[i];
            IPotential2[] p1 = atomPotentials[a1.getType().getIndex()];
            for (int j=i+1; j<atoms.size(); j++) {
                IAtom a2 = atoms.get(j);
                if (bondingInfo.skipBondedPair(false, a1, a2)) continue;
                IPotential2 p2 = p1[a2.getType().getIndex()];
                if (p2 == null) continue;
                Vector dr = space.makeVector();
                dr.Ev1Mv2(r[j], r1);
                double uu = p2.u(dr.squared());
                u += uu;
            }
        }
        return u;
    }

}
