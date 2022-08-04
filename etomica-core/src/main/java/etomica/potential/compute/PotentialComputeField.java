/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.compute;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.*;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.potential.IPotential1;
import etomica.space.Vector;
import etomica.species.SpeciesManager;
import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;

import java.util.Arrays;

public class PotentialComputeField implements PotentialCompute {

    protected final IPotential1[] potentials;
    protected final Box box;
    protected Vector[] forces, torques;
    protected double[] uAtom;
    protected final DoubleArrayList duAtom;
    protected final IntArrayList uAtomsChanged;
    protected double energyTot = Double.NaN;

    public PotentialComputeField(SpeciesManager sm, Box box) {
        this(new IPotential1[sm.getAtomTypeCount()], box);
    }

    public PotentialComputeField(IPotential1[] p1, Box box) {
        potentials = p1;
        this.box = box;
        torques = forces = new Vector[0];
        uAtom = new double[box.getLeafList().size()];
        uAtomsChanged = new IntArrayList(16);
        duAtom = new DoubleArrayList(16);

        box.getEventManager().addListener(new BoxEventListener() {
            @Override
            public void boxMoleculeAdded(BoxMoleculeEvent e) {
                int newAtoms = e.getMolecule().getType().getLeafAtomCount();
                int nowAtoms = box.getLeafList().size();
                if (nowAtoms > uAtom.length) {
                    double[] uAtomNew = new double[nowAtoms];
                    System.arraycopy(uAtom, 0, uAtomNew, 0, nowAtoms - newAtoms);
                    uAtom = uAtomNew;
                } else {
                    Arrays.fill(uAtom, nowAtoms - newAtoms, nowAtoms, 0);
                }
            }

            @Override
            public void boxNumberMolecules(BoxMoleculeCountEvent e) {
                int n = e.getCount();

                int nowAtoms = box.getLeafList().size();
                int newAtoms = e.getSpecies().getLeafAtomCount() * n;
                if (nowAtoms + newAtoms > uAtom.length) {
                    double[] uAtomNew = new double[nowAtoms + newAtoms];
                    System.arraycopy(uAtom, 0, uAtomNew, 0, nowAtoms);
                    uAtom = uAtomNew;
                }
            }

            @Override
            public void boxAtomLeafIndexChanged(BoxAtomIndexEvent e) {
                int oldIndex = e.getIndex();
                int newIndex = e.getAtom().getLeafIndex();
                uAtom[newIndex] = uAtom[oldIndex];
            }

            @Override
            public void boxMoleculeRemoved(BoxMoleculeEvent e) {
            }
        });

    }

    public void setFieldPotential(AtomType atomType, IPotential1 p) {
        potentials[atomType.getIndex()] = p;
    }

    public IPotential1[] getFieldPotentials() {
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
    public Vector[] getTorques() {
        return torques;
    }

    @Override
    public double getLastVirial() {
        return 0;
    }

    @Override
    public double getLastEnergy() {
        return energyTot;
    }

    @Override
    public void updateAtom(IAtom atom) {

    }

    protected final void zeroArrays(boolean doForces) {

        int numAtoms = box.getLeafList().size();
        if (doForces && numAtoms > forces.length) {
            int oldLength = forces.length;
            forces = Arrays.copyOf(forces, numAtoms);
            torques = Arrays.copyOf(torques, numAtoms);
            for (int i = oldLength; i < numAtoms; i++) {
                forces[i] = box.getSpace().makeVector();
                torques[i] = box.getSpace().makeVector();
            }
        }
        if (numAtoms > uAtom.length) {
            uAtom = new double[numAtoms];
        }
        for (int i = 0; i < numAtoms; i++) {
            uAtom[i] = 0;
            if (doForces) {
                forces[i].E(0);
                torques[i].E(0);
            }
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
            IPotential1 ip = potentials[iType];
            if (ip == null) continue;
            double u;
            if (doForces) {
                u = ip.uduTorque(iAtom, forces[i], torques[i]);
            } else {
                u = ip.u(iAtom);
            }
            uAtom[i] = u;
            uTot += u;
        }
        energyTot = uTot;
        return uTot;
    }

    @Override
    public double computeOneOld(IAtom iAtom) {
        return uAtom[iAtom.getLeafIndex()];
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
        IPotential1 ip = potentials[iType];
        if (ip == null) return 0;
        uAtomsChanged.add(iAtom.getLeafIndex());
        double u = ip.u(iAtom);
        duAtom.add(u);
        return u;
    }

    @Override
    public double computeManyAtomsOld(IAtom... atoms) {
        double u = 0;
        for (IAtom atom : atoms) {
            u += uAtom[atom.getLeafIndex()];
        }
        return u;
    }

    @Override
    public double computeManyAtoms(IAtom... atoms) {
        uAtomsChanged.clear();
        duAtom.clear();
        duAtom.ensureCapacity(box.getLeafList().size());
        double u = 0;

        for (IAtom atom : atoms) {
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
