/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.compute;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.*;
import etomica.integrator.IntegratorListener;
import etomica.nbr.cell.NeighborIteratorCell;
import etomica.potential.BondingInfo;
import etomica.potential.IPotential2;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space1d.Vector1D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesManager;
import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;

import java.util.Arrays;
import java.util.Objects;

public class PotentialComputePairGeneral implements PotentialCompute {

    protected final boolean isPureAtoms;
    private final NeighborIterator neighborIterator;
    protected final IPotential2[][] pairPotentials;
    protected final Box box;
    private final NeighborManager neighborManager;
    private final NeighborIteratorCell.SuperNbrConsumer nbrConsumer;
    protected double[] uAtom;
    protected final DoubleArrayList duAtom;
    protected final IntArrayList uAtomsChanged;
    protected double virialTot = Double.NaN, energyTot = Double.NaN;
    protected Vector[] forces, torques;
    protected final Space space;

    protected final int[] atomCountByType;
    protected boolean duAtomMulti = false;

    public boolean first = true;
    public boolean doAllTruncationCorrection = true;
    public boolean doOneTruncationCorrection = false;

    public PotentialComputePairGeneral(SpeciesManager sm, Box box, NeighborManager neighborManager) {
        this(sm, new IPotential2[sm.getAtomTypeCount()][sm.getAtomTypeCount()], box, neighborManager);
    }

    public PotentialComputePairGeneral(SpeciesManager sm, IPotential2[][] p2, Box box, NeighborManager neighborManager) {
        isPureAtoms = sm.isPureAtoms();
        space = box.getSpace();
        pairPotentials = p2;
        this.neighborManager = neighborManager;
        this.neighborManager.setPairPotentials(pairPotentials);
        this.neighborIterator = neighborManager.makeNeighborIterator();
        this.box = box;

        uAtom = new double[box.getLeafList().size()];
        uAtomsChanged = new IntArrayList(16);
        duAtom = new DoubleArrayList(16);
        torques = forces = new Vector[0];

        this.atomCountByType = new int[p2.length];
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

        this.nbrConsumer = new NeighborIteratorCell.SuperNbrConsumer() {
            @Override
            public double accept(IAtom atom1, IAtom jAtom, Vector rij, int n) {
                IPotential2 pij = pairPotentials[atom1.getType().getIndex()][jAtom.getType().getIndex()];
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
    public boolean needForcesForVirial() {
        return true;
    }

    @Override
    public void init() {
        this.neighborManager.init();
    }

    @Override
    public Vector[] getForces() {
        return forces;
    }

    public Vector[] getTorques() {
        return torques;
    }

    @Override
    public double getLastVirial() {
        return virialTot;
    }

    @Override
    public double getLastEnergy() {
        return energyTot;
    }

    public void setPairPotential(AtomType atomType1, AtomType atomType2, IPotential2 p12) {
        pairPotentials[atomType1.getIndex()][atomType2.getIndex()] = p12;
        pairPotentials[atomType2.getIndex()][atomType1.getIndex()] = p12;

        double maxRange = Arrays.stream(pairPotentials).flatMap(Arrays::stream)
                .filter(Objects::nonNull)
                .mapToDouble(IPotential2::getRange)
                .max().orElse(0);

        this.neighborManager.setPotentialRange(maxRange);
    }

    public IPotential2[][] getPairPotentials() {
        return pairPotentials;
    }

    @Override
    public void updateAtom(IAtom atom) {
        this.neighborManager.updateAtom(atom);
    }

    protected Vector makeTorque() {
        return box.getSpace().D() == 2 ? new Vector1D() : new Vector3D();
    }

    protected final void zeroArrays(boolean doForces) {
        virialTot = 0;

        int numAtoms = box.getLeafList().size();
        if (doForces && numAtoms > forces.length) {
            int oldLength = forces.length;
            forces = Arrays.copyOf(forces, numAtoms);
            for (int i = oldLength; i < numAtoms; i++) forces[i] = box.getSpace().makeVector();
            torques = Arrays.copyOf(torques, numAtoms);
            for (int i = oldLength; i < numAtoms; i++) torques[i] = makeTorque();
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
        double[] uAtomOld = new double[uAtom.length];
        boolean debug = false;
        if (debug) System.arraycopy(uAtom, 0, uAtomOld, 0, uAtom.length);
        zeroArrays(doForces);

        IAtomList atoms = box.getLeafList();
        double[] uTot = {0};
        for (int i = 0; i < atoms.size(); i++) {
            IAtom iAtom = atoms.get(i);
            int iType = iAtom.getType().getIndex();
            IPotential2[] ip = pairPotentials[iType];
            int finalI = i;
            neighborIterator.iterUpNeighbors(i, (jAtom, rij, n) -> {
                int j = jAtom.getLeafIndex();
                if (pc != null && pc.skipPair(finalI, j)) return;
                int jType = jAtom.getType().getIndex();
                IPotential2 pij = ip[jType];
                if (pij == null) return;
                double uij;
                if (doForces) {
                    Vector fj = space.makeVector();
                    Vector ti = makeTorque();
                    Vector tj = makeTorque();
                    fj.E(forces[j]);
                    ti.E(torques[finalI]);
                    tj.E(torques[j]);
                    uij = pij.uduTorque(rij, iAtom, jAtom, forces[finalI], forces[j], torques[finalI], torques[j]);
                    if (uij == 0) return;
                    if (pc != null) {
                        fj.ME(forces[j]);
                        ti.ME(torques[finalI]);
                        ti.TE(-1);
                        tj.ME(torques[j]);
                        tj.TE(-1);
                        pc.pairComputeGeneral(pij, iAtom, jAtom, rij, fj, ti, tj);
                    }
                    virialTot += fj.dot(rij);
                } else {
                    uij = pij.u(rij, iAtom, jAtom);
                    if (uij == 0) return;
                }
                uAtom[finalI] += 0.5 * uij;
                uAtom[j] += 0.5 * uij;
                uTot[0] += uij;
            });
        }

        double[] uCorrection = new double[1];
        double[] duCorrection = new double[1];
        computeAllTruncationCorrection(uCorrection, duCorrection);
        uTot[0] += uCorrection[0];
        virialTot += duCorrection[0];
        if (doForces && !isPureAtoms) {
            virialTot += PotentialCompute.computeVirialIntramolecular(forces, box);
        }

        if (debug && uAtom.length == uAtomOld.length && !first) {
            boolean success = true;
            for (int i = 0; i < uAtom.length; i++) {
                if (Math.abs(uAtom[i] - uAtomOld[i]) > 1e-9) {
                    if (success) System.out.println("box " + box);
                    System.out.println("uAtom diff " + i + " " + uAtom[i] + " " + uAtomOld[i]);
                    success = false;
                }
            }
            if (!success) throw new RuntimeException("oops");
        }
        first = false;
        energyTot = uTot[0];
        return uTot[0];
    }

    @Override
    public double computeOneOld(IAtom iAtom) {
        double u = this.computeOneTruncationCorrection(iAtom.getLeafIndex());
        u += uAtom[iAtom.getLeafIndex()] * 2;
        return u;
    }

    @Override
    public double computeManyAtomsOld(IAtom... atoms) {
        double u = 0;
        for (IAtom atom : atoms) {
            u += uAtom[atom.getLeafIndex()] * 2;
            u += computeOneTruncationCorrection(atom.getLeafIndex());
        }
        u -= computeIntraAtoms(atoms);
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
        double u = this.computeOneInternal(iAtom);
        u += this.computeOneTruncationCorrection(iAtom.getLeafIndex());
        return u;
    }

    protected double computeOneInternal(IAtom iAtom) {
        return this.neighborIterator.iterAndSumAllNeighbors(iAtom, this.nbrConsumer);
    }

    @Override
    public double computeManyAtoms(IAtom... atoms) {
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

        u -= computeIntraAtoms(atoms);

        return u;
    }

    protected double computeIntraAtoms(IAtom... atoms) {
        BondingInfo bondingInfo = neighborManager.getBondingInfo();
        double uIntra = 0;
        for (int i = 0; i < atoms.length; i++) {
            IAtom atom1 = atoms[i];
            for (int j = i + 1; j < atoms.length; j++) {
                IAtom atom2 = atoms[j];
                if (bondingInfo.skipBondedPair(false, atom1, atom2)) continue;
                IPotential2 pij = pairPotentials[atom1.getType().getIndex()][atom2.getType().getIndex()];
                if (pij == null) continue;
                Vector rij = space.makeVector();
                rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                box.getBoundary().nearestImage(rij);
                double uij = pij.u(rij, atom1, atom2);
                duAtom.plusEquals(atom1.getLeafIndex(), -0.5 * uij);
                duAtom.plusEquals(atom2.getLeafIndex(), -0.5 * uij);
                uIntra += uij;
            }
        }
        return uIntra;
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
                IPotential2 p = pairPotentials[i][j];
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
                p.u01TruncationCorrection(space, u, du);
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
            IPotential2 p = pairPotentials[iType][j];
            double pairDensity;
            if (iType == j) {
                pairDensity = (atomCountByType[j] - 1) / box.getBoundary().volume();
            } else {
                pairDensity = atomCountByType[j] / box.getBoundary().volume();
            }
            double integral = p.integral(space, p.getRange());
            uCorrection += pairDensity * integral;
        }
        return uCorrection;
    }
}
