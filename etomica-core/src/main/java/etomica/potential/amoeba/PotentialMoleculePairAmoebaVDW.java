/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.amoeba;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotential2;
import etomica.potential.IPotentialMolecular;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.util.collections.IntArrayList;

public class PotentialMoleculePairAmoebaVDW implements IPotentialMolecular {
    protected final Space space;
    protected final IPotential2[][] atomPotentials;
    protected final double[] reduction;
    protected final IntArrayList[][] bonding;

    public PotentialMoleculePairAmoebaVDW(Space space, SpeciesManager sm, IntArrayList[][] bonding) {
        this(space, makeAtomPotentials(sm), bonding);
    }

    public PotentialMoleculePairAmoebaVDW(Space space, IPotential2[][] atomPotentials, IntArrayList[][] bonding) {
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

    @Override
    public double energy(IMoleculeList molecules) {
        return energy(molecules.get(0), molecules.get(1));
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
//                System.out.println(a1.getLeafIndex()+" "+a2.getLeafIndex()+" "+dr+" "+Math.sqrt(dr.squared())+" "+p2.u(dr.squared()));
                u += p2.u(dr.squared());
            }
        }
        return u;
    }

}
