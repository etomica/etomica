/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.amoeba;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.potential.IPotential2;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.util.collections.IntArrayList;

public class PotentialMoleculeAmoebaVDW {
    protected final Space space;
    protected final IPotential2[][] atomPotentials;
    protected final double[] reduction;
    protected final IntArrayList[][] bonding;

    public PotentialMoleculeAmoebaVDW(Space space, SpeciesManager sm, IntArrayList[][] bonding) {
        this(space, makeAtomPotentials(sm), bonding);
    }

    public PotentialMoleculeAmoebaVDW(Space space, IPotential2[][] atomPotentials, IntArrayList[][] bonding) {
        this.space = space;
        this.atomPotentials = atomPotentials;
        reduction = new double[atomPotentials.length];
        this.bonding = bonding;
    }

    private static IPotential2[][] makeAtomPotentials(SpeciesManager sm) {
        // we could try to store the potentials more compactly, but it doesn't really matter
        ISpecies species = sm.getSpecies(sm.getSpeciesCount() - 1);
        int lastTypeIndex = species.getAtomType(species.getUniqueAtomTypeCount() - 1).getIndex();
        return new IPotential2[lastTypeIndex + 1][lastTypeIndex + 1];
    }

    public void setAtomPotential(AtomType atomType1, AtomType atomType2, IPotential2 p2) {
        atomPotentials[atomType1.getIndex()][atomType2.getIndex()] = p2;
        atomPotentials[atomType2.getIndex()][atomType1.getIndex()] = p2;
    }

    public IPotential2[][] getAtomPotentials() {
        return atomPotentials;
    }

    public void setReduction(AtomType type, double r) {
        reduction[type.getIndex()] = r;
    }

    protected Vector getReducedPosition(IAtom a) {
        Vector r1 = space.makeVector();
        r1.E(a.getPosition());
        double red = reduction[a.getType().getIndex()];
        if (red != 0) {
            IntArrayList b = bonding[a.getParentGroup().getType().getIndex()][a.getIndex()];
            if (b.size() == 1) {
                r1.TE(red);
                r1.PEa1Tv1(1-red, a.getParentGroup().getChildList().get(b.getInt(0)).getPosition());
            }
        }
        return r1;
    }

}
