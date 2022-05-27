/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.amoeba;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.integrator.IntegratorListener;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;

public class PotentialComputeIntramolecular implements PotentialCompute {

    protected final Box box;
    protected final IPotentialMoleculeSingle[] potentials;
    protected double uLast;

    public PotentialComputeIntramolecular(SpeciesManager sm, Box box) {
        potentials = new IPotentialMoleculeSingle[sm.getSpeciesCount()];
        this.box = box;
    }

    public void setPotential(ISpecies species, IPotentialMoleculeSingle p) {
        potentials[species.getIndex()] = p;
    }

    @Override
    public void init() {

    }

    @Override
    public Vector[] getForces() {
        return null;
    }

    @Override
    public double getLastVirial() {
        return 0;
    }

    @Override
    public double getLastEnergy() {
        return 0;
    }

    @Override
    public void updateAtom(IAtom atom) {

    }

    @Override
    public double computeAll(boolean doForces, PotentialCallback pc) {
        if (doForces) throw new RuntimeException("don't know how to compute forces");
        double u = 0;
        IMoleculeList molecules = box.getMoleculeList();
        for (IMolecule m : molecules) {
            u += computeOneMolecule(m);
        }
        return u;
    }

    public double computeOneOldMolecule(IMolecule molecule) {
        return computeOneMolecule(molecule);
    }

    public double computeOneMolecule(IMolecule molecule) {
        IPotentialMoleculeSingle p = potentials[molecule.getType().getIndex()];
        if (p == null) return 0;
        return p.energy(molecule);
    }


    @Override
    public double computeOneOld(IAtom iAtom) {
        throw new RuntimeException("nope");
    }

    @Override
    public double computeOne(IAtom iAtom) {
        throw new RuntimeException("nope");
    }

    @Override
    public double computeManyAtomsOld(IAtom... atoms) {
        throw new RuntimeException("nope");
    }

    @Override
    public double computeManyAtoms(IAtom... atoms) {
        throw new RuntimeException("nope");
    }

    @Override
    public void processAtomU(double fac) {

    }

    @Override
    public IntegratorListener makeIntegratorListener() {
        return new IntegratorListener() {};
    }
}
