/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesManager;

import static etomica.potential.PotentialMasterBonding.FullBondingInfo;

public class PotentialMoleculeIntra implements IPotentialMolecular {
    protected final Space space;
    protected final boolean isPureAtoms;
    protected final FullBondingInfo bondingInfo;
    protected Boundary boundary;
    protected final Potential2Soft[][] atomPotentials;

    public PotentialMoleculeIntra(Space space, SpeciesManager sm, FullBondingInfo bondingInfo) {
        this.bondingInfo = bondingInfo;
        this.space = space;
        // the species we apply to might be purely atomic.  so long as our
        // BondingInfo is nonBonding, it won't matter
        isPureAtoms = sm.isPureAtoms();
        int numAtomTypes = sm.getAtomTypeCount();
        atomPotentials = new Potential2Soft[numAtomTypes][numAtomTypes];
    }

    @Override
    public double energy(IMoleculeList molecules) {
        return energy(molecules.get(0));
    }

    public double energy(IMolecule molecule) {
        IAtomList atoms = molecule.getChildList();
        double u = 0;
        for (int i = 0; i < atoms.size(); i++) {
            IAtom a0 = atoms.get(i);
            Potential2Soft[] p0 = atomPotentials[a0.getType().getIndex()];
            for (int j = i + 1; j < atoms.size(); j++) {
                IAtom a1 = atoms.get(j);
                if (bondingInfo.skipBondedPair(isPureAtoms, a0, a1)) continue;
                Potential2Soft p2 = p0[a1.getType().getIndex()];
                if (p2 == null) continue;
                Vector dr = space.makeVector();
                dr.Ev1Mv2(a1.getPosition(), a0.getPosition());
                boundary.nearestImage(dr);
                u += p2.u(dr.squared());
            }
        }

        return u + PotentialMasterBonding.computeOneMolecule(boundary, molecule, bondingInfo);
    }

    @Override
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    @Override
    public void setBox(Box box) {
        this.boundary = box.getBoundary();
    }

    @Override
    public int nBody() {
        return 1;
    }
}
