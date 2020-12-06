/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.compute;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.molecule.IMolecule;
import etomica.potential.IPotentialAtom;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;

import java.util.Arrays;

/**
 * PotentialCompute that handles field potentials that apply to single atoms.
 * This currently does not make use of any scheme like cell-listing to speed up
 * identification of atoms that interact with the field.
 */
public class PotentialComputeOne implements PotentialCompute {

    private Vector[] forces;
    private final IPotentialAtom[] potentials;
    private final Box box;
    protected final IntArrayList uAtomsChanged;
    protected final DoubleArrayList duAtom;
    private double[] uAtom;

    PotentialComputeOne(Simulation sim, Box box) {
        ISpecies species = sim.getSpecies(sim.getSpeciesCount() - 1);
        int lastTypeIndex = species.getAtomType(species.getUniqueAtomTypeCount() - 1).getIndex();
        potentials = new IPotentialAtom[lastTypeIndex + 1];
        this.box = box;
        uAtom = new double[box.getLeafList().size()];
        uAtomsChanged = new IntArrayList(16);
        duAtom = new DoubleArrayList(16);
        forces = new Vector[0];
    }

    public void setPotential(AtomType atomType1, IPotentialAtom p1) {
        potentials[atomType1.getIndex()] = p1;
    }

    @Override
    public void init() {

    }

    @Override
    public Vector[] getForces() {
        return forces;
    }

    @Override
    public double getLastVirial() {
        return 0;
    }

    @Override
    public double getOldEnergy() {
        double uTot = 0;
        for (double iuAtom : uAtom) {
            uTot += iuAtom;
        }
        return uTot;
    }

    @Override
    public void updateAtom(IAtom atom) {

    }

    protected final void zeroArrays(boolean doForces) {
        int numAtoms = box.getLeafList().size();
        if (doForces && numAtoms > forces.length) {
            int oldLength = forces.length;
            forces = Arrays.copyOf(forces, numAtoms);
            for (int i = oldLength; i < numAtoms; i++) forces[i] = box.getSpace().makeVector();
        }
        if (numAtoms > uAtom.length) {
            uAtom = new double[numAtoms];
        }
        for (int i = 0; i < numAtoms; i++) {
            uAtom[i] = 0;
            if (doForces) forces[i].E(0);
        }
    }

    @Override
    public double computeAll(boolean doForces) {
        zeroArrays(doForces);
        IAtomList atoms = box.getLeafList();
        double uTot = 0;
        for (int i = 0; i < atoms.size(); i++) {
            IAtom iAtom = atoms.get(i);
            int iType = iAtom.getType().getIndex();
            IPotentialAtom p = potentials[iType];
            if (p == null) continue;
            double u;
            if (doForces) {
                u = p.udu(iAtom, forces[i]);
            } else {
                u = p.u(iAtom);
            }
            uAtom[i] = u;
            uTot += u;
        }

        return uTot;
    }

    @Override
    public double computeOneOld(IAtom iAtom) {
        return uAtom[iAtom.getLeafIndex()];
    }

    @Override
    public double computeOneOldMolecule(IMolecule molecule) {
        double u = 0;
        for (IAtom atom : molecule.getChildList()) {
            u += uAtom[atom.getLeafIndex()];
        }
        return u;
    }

    @Override
    public double computeOne(IAtom iAtom) {
        int i = iAtom.getLeafIndex();
        IAtomList atoms = box.getLeafList();
        uAtomsChanged.clear();
        uAtomsChanged.ensureCapacity(atoms.size());
        uAtomsChanged.add(i);
        duAtom.clear();
        duAtom.ensureCapacity(atoms.size());
        duAtom.add(0);
        return this.computeOneInternal(iAtom);
    }

    private double computeOneInternal(IAtom iAtom) {
        int iType = iAtom.getType().getIndex();
        IPotentialAtom p = potentials[iType];
        if (p == null) return 0;
        double u = p.u(iAtom);
        uAtomsChanged.add(iAtom.getLeafIndex());
        duAtom.add(u);
        return u;
    }

    @Override
    public double computeOneMolecule(IMolecule molecule) {
        double u = 0;
        for (IAtom atom : molecule.getChildList()) {
            u += computeOneInternal(atom);
        }
        return u;
    }

    @Override
    public void processAtomU(double fac) {
        for (int j = 0; j < uAtomsChanged.size(); j++) {
            int jj = uAtomsChanged.getInt(j);
            uAtom[jj] += fac * duAtom.getDouble(j);
        }
    }

    @Override
    public IntegratorListener makeIntegratorListener() {
        return new IntegratorListener() {
            @Override
            public void integratorInitialized(IntegratorEvent e) {

            }

            @Override
            public void integratorStepStarted(IntegratorEvent e) {

            }

            @Override
            public void integratorStepFinished(IntegratorEvent e) {

            }
        };
    }
}
