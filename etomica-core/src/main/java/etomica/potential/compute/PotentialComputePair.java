/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.compute;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.*;
import etomica.integrator.IntegratorListener;
import etomica.nbr.cell.NeighborIteratorCellFasterer;
import etomica.potential.BondingInfo;
import etomica.potential.IPotential;
import etomica.potential.Potential2Soft;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesManager;
import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;

import java.util.Arrays;
import java.util.Objects;

public class PotentialComputePair implements PotentialCompute {
    private final NeighborIterator neighborIterator;
    protected final Potential2Soft[][] pairPotentials;
    protected final Box box;
    private final NeighborManager neighborManager;
    private final NeighborIteratorCellFasterer.SuperNbrConsumer nbrConsumer;
    protected double[] uAtom;
    protected final DoubleArrayList duAtom;
    protected final IntArrayList uAtomsChanged;
    protected Vector zero;
    protected double virialTot = Double.NaN, energyTot = Double.NaN;
    protected Vector[] forces;
    protected final Space space;

    protected final int[] atomCountByType;
    protected boolean duAtomMulti = false;

    public boolean doAllTruncationCorrection = true;
    public boolean doOneTruncationCorrection = false;

    public PotentialComputePair(SpeciesManager sm, Box box, NeighborManager neighborManager) {
        this(sm, box, neighborManager, new Potential2Soft[sm.getAtomTypeCount()][sm.getAtomTypeCount()]);
    }

    public PotentialComputePair(SpeciesManager sm, Box box, NeighborManager neighborManager, Potential2Soft[][] pairPotentials) {
        space = box.getSpace();
        int numAtomTypes = sm.getAtomTypeCount();
        this.neighborManager = neighborManager;
        this.pairPotentials = pairPotentials;
        updateNeighborRange();

        this.neighborManager.setPairPotentials(pairPotentials);
        this.neighborIterator = neighborManager.makeNeighborIterator();
        this.box = box;

        uAtom = new double[box.getLeafList().size()];
        uAtomsChanged = new IntArrayList(16);
        duAtom = new DoubleArrayList(16);
        zero = box.getSpace().makeVector();
        forces = new Vector[0];

        this.atomCountByType = new int[numAtomTypes];
        box.getEventManager().addListener(new BoxEventListener() {
            @Override
            public void boxMoleculeAdded(BoxMoleculeEvent e) {
                for (AtomType atomType : e.getMolecule().getType().getAtomTypes()) {
                    atomCountByType[atomType.getIndex()]++;
                }

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
                for (AtomType atomType : e.getMolecule().getType().getAtomTypes()) {
                    atomCountByType[atomType.getIndex()]--;
                }
            }
        });

        this.nbrConsumer = new NeighborIteratorCellFasterer.SuperNbrConsumer() {
            @Override
            public double accept(IAtom atom1, IAtom jAtom, Vector rij) {
                Potential2Soft pij = pairPotentials[atom1.getType().getIndex()][jAtom.getType().getIndex()];
                if (pij == null) return 0;
                double uij = pij.u(rij.squared());
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

    public Potential2Soft[][] getPairPotentials() {
        return pairPotentials;
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

    public void setPairPotential(AtomType atomType1, AtomType atomType2, Potential2Soft p12) {
        pairPotentials[atomType1.getIndex()][atomType2.getIndex()] = p12;
        pairPotentials[atomType2.getIndex()][atomType1.getIndex()] = p12;
        updateNeighborRange();
    }

    protected void updateNeighborRange() {
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
            Potential2Soft[] ip = pairPotentials[iType];
            int finalI = i;
            neighborIterator.iterUpNeighbors(i, (jAtom, rij) -> {
                int j = jAtom.getLeafIndex();
                if (pc != null && pc.skipPair(finalI, j)) return;
                int jType = jAtom.getType().getIndex();
                Potential2Soft pij = ip[jType];
                double[] u012 = new double[3];
                double r2 = rij.squared();
                pij.u012add(r2, u012);
                double uij = u012[0];
                if (uij == 0) return;
                if (pc != null) pc.pairCompute(finalI, j, rij, u012);
                uAtom[finalI] += 0.5 * uij;
                uAtom[j] += 0.5 * uij;
                double duij = u012[1];
                virialTot += duij;
                if (doForces) {
                    rij.TE(duij / r2);
                    forces[finalI].PE(rij);
                    forces[j].ME(rij);
                }
                uTot[0] += uij;
            });
        }

        double[] uCorrection = new double[1];
        double[] duCorrection = new double[1];
        this.computeAllTruncationCorrection(uCorrection, duCorrection);
        uTot[0] += uCorrection[0];
        virialTot += duCorrection[0];
        energyTot = uTot[0];
        return uTot[0];
    }

    @Override
    public double computeOneOld(IAtom iAtom) {
        double u = this.computeOneTruncationCorrection(iAtom.getLeafIndex());
        return u + uAtom[iAtom.getLeafIndex()] * 2;
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
        double u = this.computeOneInternal(iAtom);
        u += this.computeOneTruncationCorrection(i);
        return u;
    }

    protected double computeOneInternal(IAtom iAtom) {
        return this.neighborIterator.iterAndSumAllNeighbors(iAtom, this.nbrConsumer);
    }

    @Override
    public double computeManyAtomsOld(IAtom... atoms) {
        double u = 0;
        for (IAtom atom : atoms) {
            u += uAtom[atom.getLeafIndex()] * 2;
            u += this.computeOneTruncationCorrection(atom.getLeafIndex());
        }

        u -= computeIntraAtoms(true, atoms);

        return u;
    }

    @Override
    public double computeManyAtoms(IAtom... atoms) {
        return computeManyAtoms(true, atoms);
    }

    protected double computeIntraAtoms(boolean doOld, IAtom... atoms) {
        BondingInfo bondingInfo = neighborManager.getBondingInfo();
        double uIntra = 0;
        for (int i = 0; i < atoms.length; i++) {
            IAtom atom1 = atoms[i];
            for (int j = i + 1; j < atoms.length; j++) {
                IAtom atom2 = atoms[j];
                if (bondingInfo.skipBondedPair(false, atom1, atom2)) continue;
                Potential2Soft pij = pairPotentials[atom1.getType().getIndex()][atom2.getType().getIndex()];
                if (pij == null) continue;
                Vector rij = space.makeVector();
                rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                box.getBoundary().nearestImage(rij);
                double uij = pij.u(rij.squared());
                if (!doOld) {
                    duAtom.plusEquals(atom1.getLeafIndex(), -0.5 * uij);
                    duAtom.plusEquals(atom2.getLeafIndex(), -0.5 * uij);
                }
                uIntra += uij;
            }
        }
        return uIntra;
    }

    protected double computeManyAtoms(boolean doInternal, IAtom... atoms) {
        duAtomMulti = true;
        uAtomsChanged.clear();
        duAtom.setAll(0);
        duAtom.ensureCapacity(box.getLeafList().size());
        double u = 0;

        for (IAtom atom : atoms) {
            uAtomsChanged.add(atom.getLeafIndex());

            u += computeOneInternal(atom);
            u += computeOneTruncationCorrection(atom.getLeafIndex());
        }

        if (doInternal) {
            u -= computeIntraAtoms(false, atoms);
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

    public void computeAllTruncationCorrection(double[] uCorrection, double[] duCorrection) {
        if (!doAllTruncationCorrection) {
            return;
        }
        double uCor = 0;
        double duCor = 0;
        for (int i = 0; i < atomCountByType.length; i++) {
            for (int j = i; j < atomCountByType.length; j++) {
                Potential2Soft p = pairPotentials[i][j];
                if (p == null) continue;
                int numPairs;
                if (j == i) {
                    numPairs = atomCountByType[i] * (atomCountByType[j] - 1) / 2;
                } else {
                    numPairs = atomCountByType[i] * atomCountByType[j];
                }
                double pairDensity = numPairs / box.getBoundary().volume();

                double[] u = new double[1];
                double[] du = new double[1];
                p.u01TruncationCorrection(u, du);
                uCor += pairDensity * u[0];
                duCor += pairDensity * du[0];

            }
        }
        uCorrection[0] = uCor;
        duCorrection[0] = duCor;
    }

    public double computeOneTruncationCorrection(int iAtom) {
        if (!doOneTruncationCorrection) {
            return 0;
        }
        int iType = box.getLeafList().get(iAtom).getType().getIndex();
        double uCorrection = 0;
        for (int j = 0; j < atomCountByType.length; j++) {
            Potential2Soft p = pairPotentials[iType][j];
            double pairDensity;
            if (iType == j) {
                pairDensity = (atomCountByType[j] - 1) / box.getBoundary().volume();
            } else {
                pairDensity = atomCountByType[j] / box.getBoundary().volume();
            }
            double integral = p.integral(p.getRange());
            uCorrection += pairDensity * integral;
        }
        return uCorrection;
    }
}
