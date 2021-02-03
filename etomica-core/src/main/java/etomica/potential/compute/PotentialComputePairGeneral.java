/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.compute;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.box.BoxEventListener;
import etomica.box.BoxMoleculeEvent;
import etomica.integrator.IntegratorListener;
import etomica.molecule.IMolecule;
import etomica.nbr.cell.NeighborIteratorCellFasterer;
import etomica.potential.IPotential;
import etomica.potential.IPotentialPair;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;

import java.util.Arrays;
import java.util.Objects;

public class PotentialComputePairGeneral implements PotentialCompute {

    private final NeighborIterator neighborIterator;
    protected final IPotentialPair[][] pairPotentials;
    protected final Box box;
    private final NeighborManager neighborManager;
    private final NeighborIteratorCellFasterer.SuperNbrConsumer nbrConsumer;
    protected double[] uAtom;
    protected final DoubleArrayList duAtom;
    protected final IntArrayList uAtomsChanged;
    protected double virialTot = Double.NaN, energyTot = Double.NaN;
    protected Vector[] forces;
    protected final Space space;

    protected final int[] atomCountByType;
    protected boolean duAtomMulti = false;

    public PotentialComputePairGeneral(Simulation sim, Box box, NeighborManager neighborManager) {
        space = box.getSpace();
        ISpecies species = sim.getSpecies(sim.getSpeciesCount() - 1);
        int lastTypeIndex = species.getAtomType(species.getUniqueAtomTypeCount() - 1).getIndex();
        pairPotentials = new IPotentialPair[lastTypeIndex + 1][lastTypeIndex + 1];
        this.neighborManager = neighborManager;
        this.neighborManager.setPairPotentials(pairPotentials);
        this.neighborIterator = neighborManager.makeNeighborIterator();
        this.box = box;

        uAtom = new double[box.getLeafList().size()];
        uAtomsChanged = new IntArrayList(16);
        duAtom = new DoubleArrayList(16);
        forces = new Vector[0];

        this.atomCountByType = new int[lastTypeIndex + 1];
        box.getEventManager().addListener(new BoxEventListener() {
            @Override
            public void boxMoleculeAdded(BoxMoleculeEvent e) {
                for (AtomType atomType : e.getMolecule().getType().getAtomTypes()) {
                    atomCountByType[atomType.getIndex()]++;
                }
            }

            @Override
            public void boxMoleculeRemoved(BoxMoleculeEvent e) {
                for (AtomType atomType : e.getMolecule().getType().getAtomTypes()) {
                    atomCountByType[atomType.getIndex()]--;
                }
            }
        });

        this.nbrConsumer = new NeighborIteratorCellFasterer.SuperNbrConsumer() {
            @Override
            public double accept(IAtom atom1, IAtom jAtom, Vector rij) {
                IPotentialPair pij = pairPotentials[atom1.getType().getIndex()][jAtom.getType().getIndex()];
                if (pij == null) return 0;
                double uij = pij.u(rij, atom1, jAtom);
                uAtomsChanged.add(jAtom.getLeafIndex());
                if (duAtomMulti) {
                    duAtom.plusEquals(atom1.getLeafIndex(), 0.5 * uij);
                    duAtom.plusEquals(jAtom.getLeafIndex(), 0.5 * uij);
                } else {
                    duAtom.plusEquals(0, 0.5 * uij);
                    duAtom.add(0.5 * uij);
                }
                return uij;
            }
        };
    }

    @Override
    public void init() {
        this.neighborManager.init();
    }

    @Override
    public Vector[] getForces() {
        return forces;
    }

    @Override
    public double getLastVirial() {
        return virialTot;
    }

    @Override
    public double getLastEnergy() {
        return energyTot;
    }

    public void setPairPotential(AtomType atomType1, AtomType atomType2, IPotentialPair p12) {
        pairPotentials[atomType1.getIndex()][atomType2.getIndex()] = p12;
        pairPotentials[atomType2.getIndex()][atomType1.getIndex()] = p12;

        double maxRange = Arrays.stream(pairPotentials).flatMap(Arrays::stream)
                .filter(Objects::nonNull)
                .mapToDouble(IPotential::getRange)
                .max().orElse(0);

        this.neighborManager.setPotentialRange(maxRange);
    }

    @Override
    public void updateAtom(IAtom atom) {
        this.neighborManager.updateAtom(atom);
    }

    protected final void zeroArrays(boolean doForces) {
        virialTot = 0;

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
        double[] uTot = {0};
        for (int i = 0; i < atoms.size(); i++) {
            IAtom iAtom = atoms.get(i);
            int iType = iAtom.getType().getIndex();
            IPotentialPair[] ip = pairPotentials[iType];
            int finalI = i;
            neighborIterator.iterUpNeighbors(i, (jAtom, rij) -> {
                int j = jAtom.getLeafIndex();
                int jType = jAtom.getType().getIndex();
                IPotentialPair pij = ip[jType];
                double uij;
                if (doForces) {
                    uij = pij.udu(rij, iAtom, jAtom, forces[finalI], forces[j]);
                    if (uij == 0) return;
                    virialTot -= forces[j].dot(rij);
                } else {
                    uij = pij.u(rij, iAtom, jAtom);
                    if (uij == 0) return;
                }
                uAtom[finalI] += 0.5 * uij;
                uAtom[j] += 0.5 * uij;
            });
        }
        energyTot = uTot[0];
        return uTot[0];
    }

    @Override
    public double computeOneOld(IAtom iAtom) {
        return uAtom[iAtom.getLeafIndex()] * 2;
    }

    @Override
    public double computeOneOldMolecule(IMolecule molecule) {
        double u = 0;
        for (IAtom atom : molecule.getChildList()) {
            u += uAtom[atom.getLeafIndex()] * 2;
        }
        return u;
    }

    @Override
    public double computeOne(IAtom iAtom) {
        this.duAtomMulti = false;
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

    protected double computeOneInternal(IAtom iAtom) {
        return this.neighborIterator.iterAndSumAllNeighbors(iAtom, this.nbrConsumer);
    }

    @Override
    public double computeOneMolecule(IMolecule molecule) {
        duAtomMulti = true;
        uAtomsChanged.clear();
        duAtom.setAll(0);
        duAtom.ensureCapacity(box.getLeafList().size());
        double u = 0;

        for (IAtom atom : molecule.getChildList()) {
            uAtomsChanged.add(atom.getLeafIndex());

            u += computeOneInternal(atom);
        }

        return u;
    }

    @Override
    public void processAtomU(double fac) {
        for (int j = 0; j < uAtomsChanged.size(); j++) {
            int jj = uAtomsChanged.getInt(j);
            if (duAtomMulti) {
                uAtom[jj] += fac * duAtom.getDouble(jj);
                duAtom.replace(jj, 0);
            } else {
                uAtom[jj] += fac * duAtom.getDouble(j);
            }
        }
    }

    @Override
    public IntegratorListener makeIntegratorListener() {
        return this.neighborManager.makeIntegratorListener();
    }
}
