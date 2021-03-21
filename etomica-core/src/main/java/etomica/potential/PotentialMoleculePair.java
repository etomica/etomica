/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;

public class PotentialMoleculePair implements IPotentialMolecular {
    protected final Space space;
    protected Boundary boundary;
    protected final Potential2Soft[][] atomPotentials;

    public PotentialMoleculePair(Space space, SpeciesManager sm) {
        this.space = space;
        // we could try to store the potentials more compactly, but it doesn't really matter
        ISpecies species = sm.getSpecies(sm.getSpeciesCount() - 1);
        int lastTypeIndex = species.getAtomType(species.getUniqueAtomTypeCount() - 1).getIndex();
        atomPotentials = new Potential2Soft[lastTypeIndex + 1][lastTypeIndex + 1];
    }

    public void setAtomPotential(AtomType atomType1, AtomType atomType2, Potential2Soft p2) {
        atomPotentials[atomType1.getIndex()][atomType2.getIndex()] = p2;
        atomPotentials[atomType2.getIndex()][atomType1.getIndex()] = p2;
    }

    @Override
    public double energy(IMoleculeList molecules) {
        return energy(molecules.get(0), molecules.get(1));
    }

    public double energy(IMolecule molecule1, IMolecule molecule2) {
        IAtomList atoms1 = molecule1.getChildList();
        IAtomList atoms2 = molecule2.getChildList();
        double u = 0;
        for (IAtom a0 : atoms1) {
            Potential2Soft[] p0 = atomPotentials[a0.getType().getIndex()];
            for (IAtom a1 : atoms2) {
                Potential2Soft p2 = p0[a1.getType().getIndex()];
                if (p2 == null) continue;
                Vector dr = space.makeVector();
                dr.Ev1Mv2(a1.getPosition(), a0.getPosition());
                boundary.nearestImage(dr);
                u += p2.u(dr.squared());
            }
        }
        return u;
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
        return 2;
    }
}
