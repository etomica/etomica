/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;

public class PotentialMoleculePair implements IPotentialMolecular {
    protected final Space space;
    protected final IPotential2[][] atomPotentials;
    protected double[][][] pScale;

    public PotentialMoleculePair(Space space, SpeciesManager sm) {
        this(space, makeAtomPotentials(sm));
        pScale = new double[sm.getAtomTypeCount()][sm.getAtomTypeCount()][];
       // this(space, makeAtomPotentials(sm), pScale);
    }

    public PotentialMoleculePair(Space space, IPotential2[][] atomPotentials) {
        this.space = space;
        this.atomPotentials = atomPotentials;
    }
    public PotentialMoleculePair(Space space, IPotential2[][] atomPotentials, double[][][] pScale) {
        this.space = space;
        this.atomPotentials = atomPotentials;
        this.pScale = pScale;
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
    public void setAtomPotentials(AtomType atomType1, AtomType atomType2, IPotential2 p12, double[] ijScale){
        atomPotentials[atomType1.getIndex()][atomType2.getIndex()] = p12;
        atomPotentials[atomType2.getIndex()][atomType1.getIndex()] = p12;
        pScale[atomType2.getIndex()][atomType1.getIndex()]= ijScale;
        pScale[atomType1.getIndex()][atomType2.getIndex()]= ijScale;
    }

    public IPotential2[][] getAtomPotentials() {
        return atomPotentials;
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
            IPotential2[] p0 = atomPotentials[a0.getType().getIndex()];
            for (IAtom a1 : atoms2) {
                IPotential2 p2 = p0[a1.getType().getIndex()];
                if (p2 == null) continue;
                Vector dr = space.makeVector();
                dr.Ev1Mv2(a1.getPosition(), a0.getPosition());
                u += p2.u(dr.squared());
                if(u>8E7){
                  //  System.out.println("Next One");
                }
                //System.out.println(" sum : "+ u + " " + a0 +" "+a1);
            }
        }
        return u;
    }

}
