/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.*;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.potential.compute.PotentialCallback;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;

import java.util.Arrays;

public class PotentialMasterFasterer implements etomica.potential.compute.PotentialCompute {
    protected final BondingInfo bondingInfo;
    protected final Potential2Soft[][] pairPotentials;
    protected final Box box;
    protected final Vector dr;
    protected double[] uAtom;
    protected final DoubleArrayList duAtom;
    protected final IntArrayList uAtomsChanged;
    protected Vector zero;
    protected double virialTot = Double.NaN, energyTot = Double.NaN;
    protected Vector[] forces;
    protected final Space space;
    protected double minR2;

    protected final boolean isPureAtoms;
    protected final int[] atomCountByType;
    protected boolean duAtomMulti = false;

    public boolean doAllTruncationCorrection = true;
    public boolean doOneTruncationCorrection = false;
    boolean first = true;
    public long tAll, numAll;
    public long numMC, tMC;

    public PotentialMasterFasterer(SpeciesManager sm, Box box, BondingInfo bondingInfo) {
        this(sm, box, bondingInfo, new Potential2Soft[sm.getAtomTypeCount()][sm.getAtomTypeCount()]);
    }

    public PotentialMasterFasterer(SpeciesManager sm, Box box, BondingInfo bondingInfo, Potential2Soft[][] pairPotentials) {
        space = box.getSpace();
        this.bondingInfo = bondingInfo;
        this.pairPotentials = pairPotentials;
        this.box = box;
        dr = box.getSpace().makeVector();

        uAtom = new double[box.getLeafList().size()];
        uAtomsChanged = new IntArrayList(16);
        duAtom = new DoubleArrayList(16);
        zero = box.getSpace().makeVector();
        forces = new Vector[0];

        isPureAtoms = sm.isPureAtoms();

        this.atomCountByType = new int[sm.getAtomTypeCount()];
        for (ISpecies s : sm.getSpeciesArray()) {
            int nMols = box.getNMolecules(s);
            for (AtomType type : s.getAtomTypes()) {
                atomCountByType[type.getIndex()] += nMols;
            }
        }
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
        return virialTot;
    }

    @Override
    public double getLastEnergy() {
        return energyTot;
    }

    public void setPairPotential(AtomType atomType1, AtomType atomType2, Potential2Soft p12) {
        pairPotentials[atomType1.getIndex()][atomType2.getIndex()] = p12;
        pairPotentials[atomType2.getIndex()][atomType1.getIndex()] = p12;
    }

    public Potential2Soft[][] getPairPotentials() {
        return pairPotentials;
    }

    @Override
    public void updateAtom(IAtom atom) {

    }

    protected double handleComputeAll(boolean doForces, int iAtom, int jAtom, Vector ri, Vector rj, Vector jbo, Potential2Soft pij, PotentialCallback pc, boolean skipIntra) {
        if (pc != null && pc.skipPair(iAtom, jAtom)) return 0;
        numAll++;
        dr.Ev1Mv2(rj, ri);
        dr.PE(jbo);
        double[] u012 = new double[3];
        double r2 = dr.squared();
        if (skipIntra && r2 < minR2) return 0;
        pij.u012add(r2, u012);
        double uij = u012[0];
//        double uij = pij.u(dr.squared());
        if (uij == 0) return 0;
        if (pc != null) pc.pairCompute(iAtom, jAtom, dr, u012);
//        System.out.println(iAtom+" "+jAtom+" "+uij);
        uAtom[iAtom] += 0.5 * uij;
        uAtom[jAtom] += 0.5 * uij;
        double duij = u012[1];
        virialTot += duij;
        if (doForces) {
            dr.TE(duij / r2);
            forces[iAtom].PE(dr);
            forces[jAtom].ME(dr);
        }
        return uij;
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
        double[] uAtomOld = new double[uAtom.length];
        boolean debug = false;
        if (debug) System.arraycopy(uAtom, 0, uAtomOld, 0, uAtom.length);
        zeroArrays(doForces);

        IAtomList atoms = box.getLeafList();
        double u = 0;
        Boundary boundary = box.getBoundary();
        long t1 = System.nanoTime();
        for (int i = 0; i < atoms.size(); i++) {
            IAtom iAtom = atoms.get(i);
            Vector ri = iAtom.getPosition();
            int iType = iAtom.getType().getIndex();
            Potential2Soft[] ip = pairPotentials[iType];
            for (int j = i + 1; j < atoms.size(); j++) {
                IAtom jAtom = atoms.get(j);
                int jType = jAtom.getType().getIndex();
                Potential2Soft pij = ip[jType];
                if (pij == null) continue;

                if (bondingInfo.skipBondedPair(isPureAtoms, iAtom, jAtom)) continue;

                dr.Ev1Mv2(jAtom.getPosition(), ri);
                boundary.nearestImage(dr);
                u += handleComputeAll(doForces, i, j, zero, dr, zero, pij, pc, false);
            }
        }
        tAll += System.nanoTime() - t1;

        double[] uCorrection = new double[1];
        double[] duCorrection = new double[1];
        this.computeAllTruncationCorrection(uCorrection, duCorrection);
        u += uCorrection[0];
        virialTot += duCorrection[0];

        if (debug && uAtom.length == uAtomOld.length && !first) {
            boolean success = true;
            for (int i = 0; i < uAtom.length; i++) {
                if (Math.abs(uAtom[i] - uAtomOld[i]) > 1e-9) {
                    System.out.println("uAtom diff " + i + " " + uAtom[i] + " " + uAtomOld[i]);
                    success = false;
                }
            }
            if (!success) throw new RuntimeException("oops");
        }
        first = false;
        energyTot = u;
        return u;
    }

    @Override
    public double computeOneOld(IAtom iAtom) {
        double u = this.computeOneTruncationCorrection(iAtom.getLeafIndex());
        return u + uAtom[iAtom.getLeafIndex()] * 2;
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

    protected double computeIntraAtoms(boolean isOld, IAtom... atoms) {
        double uIntra = 0;
        for (int i = 0; i < atoms.length; i++) {
            IAtom atom1 = atoms[i];
            for (int j = i + 1; j < atoms.length; j++) {
                IAtom atom2 = atoms[j];
                if (bondingInfo.skipBondedPair(isPureAtoms, atom1, atom2)) continue;
                Potential2Soft pij = pairPotentials[atom1.getType().getIndex()][atom2.getType().getIndex()];
                if (pij == null) continue;
                Vector rij = space.makeVector();
                rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                box.getBoundary().nearestImage(rij);
                double uij = pij.u(rij.squared());
                if (!isOld) {
                    duAtom.plusEquals(atom1.getLeafIndex(), -0.5 * uij);
                    duAtom.plusEquals(atom2.getLeafIndex(), -0.5 * uij);
                }

                uIntra += uij;
            }
        }
        return uIntra;
    }

    protected double handleComputeOne(Potential2Soft pij, Vector ri, Vector rj, Vector jbo, int iAtom, int jAtom, boolean skipIntra) {
        numMC++;
        dr.Ev1Mv2(rj, ri);
        dr.PE(jbo);
        double r2 = dr.squared();
        if (skipIntra && r2 < minR2) return 0;
        double uij = pij.u(r2);
        if (uij == 0) return 0;

        uAtomsChanged.add(jAtom);
        if (duAtomMulti) {
            duAtom.plusEquals(iAtom, 0.5 * uij);
            duAtom.plusEquals(jAtom, 0.5 * uij);
        } else {
            duAtom.plusEquals(0, 0.5 * uij);
            duAtom.add(0.5 * uij);
        }
        return uij;
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

        double u = computeOneInternal(iAtom);
        u += this.computeOneTruncationCorrection(i);
        return u;
    }

    protected double computeOneInternal(IAtom atom) {
        int iType = atom.getType().getIndex();
        int i = atom.getLeafIndex();
        Potential2Soft[] ip = pairPotentials[iType];
        double u = 0;
        Boundary boundary = box.getBoundary();
        long t1 = System.nanoTime();
        for (int j = 0; j < box.getLeafList().size(); j++) {
            if (i == j) continue;
            IAtom jAtom = box.getLeafList().get(j);
            int jType = jAtom.getType().getIndex();
            Potential2Soft pij = ip[jType];
            if (pij == null) continue;
            if (bondingInfo.skipBondedPair(isPureAtoms, atom, jAtom)) continue;
            dr.Ev1Mv2(jAtom.getPosition(), atom.getPosition());
            boundary.nearestImage(dr);
            u += handleComputeOne(pij, zero, dr, zero, i, j, false);
        }
        tMC += System.nanoTime() - t1;
        return u;
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
        if (u == Double.POSITIVE_INFINITY) return u;

        u -= computeIntraAtoms(false, atoms);

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
        return new IntegratorListener() {
            @Override
            public void integratorInitialized(IntegratorEvent e) {
                init();
            }

            @Override
            public void integratorStepStarted(IntegratorEvent e) {

            }

            @Override
            public void integratorStepFinished(IntegratorEvent e) {

            }
        };
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
