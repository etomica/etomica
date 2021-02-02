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
import etomica.potential.IPotentialField;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;

import java.util.Arrays;

public class PotentialComputeField implements PotentialCompute {

    protected final IPotentialField[] potentials;
    protected final Box box;
    protected Vector[] forces;
    protected double[] uAtom;
    protected final DoubleArrayList duAtom;
    protected final IntArrayList uAtomsChanged;

    public PotentialComputeField(Simulation sim, Box box) {
        ISpecies species = sim.getSpecies(sim.getSpeciesCount() - 1);
        int lastTypeIndex = species.getAtomType(species.getUniqueAtomTypeCount() - 1).getIndex();
        potentials = new IPotentialField[lastTypeIndex + 1];
        this.box = box;
        forces = new Vector[0];
        uAtom = new double[box.getLeafList().size()];
        uAtomsChanged = new IntArrayList(16);
        duAtom = new DoubleArrayList(16);
    }

    public void setFieldPotential(AtomType atomType, IPotentialField p) {
        potentials[atomType.getIndex()] = p;
    }

    public IPotentialField[] getFieldPotentials() {
        return potentials;
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
    public double computeAll(boolean doForces, PotentialCallback pc) {
        zeroArrays(doForces);
        IAtomList atoms = box.getLeafList();
        double uTot = 0;
        for (int i = 0; i < atoms.size(); i++) {
            IAtom iAtom = atoms.get(i);
            int iType = iAtom.getType().getIndex();
            IPotentialField ip = potentials[iType];
            if (ip == null) continue;
            double u;
            if (doForces) {
                u = ip.udu(iAtom, forces[i]);
            } else {
                u = ip.u(iAtom);
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
        IAtomList atoms = box.getLeafList();
        uAtomsChanged.clear();
        uAtomsChanged.ensureCapacity(atoms.size());
        duAtom.clear();
        duAtom.ensureCapacity(atoms.size());
        return this.computeOneInternal(iAtom);
    }

    protected double computeOneInternal(IAtom iAtom) {
        int iType = iAtom.getType().getIndex();
        IPotentialField ip = potentials[iType];
        if (ip == null) return 0;
        uAtomsChanged.add(iAtom.getLeafIndex());
        double u = ip.u(iAtom);
        duAtom.add(u);
        return u;
    }

    @Override
    public double computeOneMolecule(IMolecule molecule) {
        uAtomsChanged.clear();
        duAtom.setAll(0);
        duAtom.ensureCapacity(box.getLeafList().size());
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
