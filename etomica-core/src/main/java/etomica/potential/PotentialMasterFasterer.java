/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.box.BoxEventListener;
import etomica.box.BoxMoleculeEvent;
import etomica.molecule.IMolecule;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;

import java.util.*;

public class PotentialMasterFasterer implements etomica.potential.compute.PotentialCompute {
    private final List<ISpecies> speciesList;
    protected final BondingInfo bondingInfo;
    protected final Potential2Soft[][] pairPotentials;
    protected final Box box;
    protected final Vector dr;
    protected double[] uAtom;
    protected final DoubleArrayList duAtom;
    protected final IntArrayList uAtomsChanged;
    protected Vector zero;
    protected double virialTot;
    protected Vector[] forces;
    protected final Space space;

    protected final boolean isPureAtoms;
    protected final int[] atomCountByType;
    protected boolean duAtomMulti = false;

    public boolean doAllTruncationCorrection = true;
    public boolean doOneTruncationCorrection = false;

    public PotentialMasterFasterer(Simulation sim, Box box, BondingInfo bondingInfo) {
        space = box.getSpace();
        this.speciesList = sim.getSpeciesList();
        this.bondingInfo = bondingInfo;
        ISpecies species = sim.getSpecies(sim.getSpeciesCount() - 1);
        int lastTypeIndex = species.getAtomType(species.getUniqueAtomTypeCount() - 1).getIndex();
        pairPotentials = new Potential2Soft[lastTypeIndex + 1][lastTypeIndex + 1];
        this.box = box;
        dr = sim.getSpace().makeVector();

        uAtom = new double[box.getLeafList().size()];
        uAtomsChanged = new IntArrayList(16);
        duAtom = new DoubleArrayList(16);
        zero = box.getSpace().makeVector();
        forces = new Vector[0];

        isPureAtoms = speciesList.stream().allMatch(s -> s.getLeafAtomCount() == 1);

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
    public double getOldEnergy() {
        double uTot = 0;
        for (double iuAtom : uAtom) {
            uTot += iuAtom;
        }
        double[] uCorrection = new double[1];
        double[] duCorrection = new double[1];
        this.computeAllTruncationCorrection(uCorrection, duCorrection);
        uTot += uCorrection[0];
        return uTot;
    }

    public void setPairPotential(AtomType atomType1, AtomType atomType2, Potential2Soft p12) {
        pairPotentials[atomType1.getIndex()][atomType2.getIndex()] = p12;
        pairPotentials[atomType2.getIndex()][atomType1.getIndex()] = p12;
    }

    @Override
    public void updateAtom(IAtom atom) {

    }

    protected double handleComputeAll(boolean doForces, int iAtom, int jAtom, Vector ri, Vector rj, Vector jbo, Potential2Soft pij) {
        dr.Ev1Mv2(rj, ri);
        dr.PE(jbo);
        double[] u = {0};
        double[] du = {0};
        double r2 = dr.squared();
        pij.udu(r2, u, du);
        double uij = u[0];
//        double uij = pij.u(dr.squared());
        if (uij == 0) return 0;
//        System.out.println(iAtom+" "+jAtom+" "+uij);
        uAtom[iAtom] += 0.5 * uij;
        uAtom[jAtom] += 0.5 * uij;
        double duij = du[0];
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
    public double computeAll(boolean doForces) {
        zeroArrays(doForces);

        IAtomList atoms = box.getLeafList();
        double u = 0;
        Boundary boundary = box.getBoundary();
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
                u += handleComputeAll(doForces, i, j, zero, dr, zero, pij);
            }
        }

        double[] uCorrection = new double[1];
        double[] duCorrection = new double[1];
        this.computeAllTruncationCorrection(uCorrection, duCorrection);
        u += uCorrection[0];
        virialTot += duCorrection[0];

        return u;
    }

    @Override
    public double computeOneOld(IAtom iAtom) {
        double u = this.computeOneTruncationCorrection(iAtom.getLeafIndex());
        return u + uAtom[iAtom.getLeafIndex()] * 2;
    }

    @Override
    public double computeOneOldMolecule(IMolecule molecule) {
        double u = 0;
        for (IAtom atom : molecule.getChildList()) {
            u += uAtom[atom.getLeafIndex()] * 2;
            u += this.computeOneTruncationCorrection(atom.getLeafIndex());
        }
        return u;
    }

    protected double handleComputeOne(Potential2Soft pij, Vector ri, Vector rj, Vector jbo, int iAtom, int jAtom) {
        dr.Ev1Mv2(rj, ri);
        dr.PE(jbo);
        double uij = pij.u(dr.squared());
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
        for (int j = 0; j < box.getLeafList().size(); j++) {
            if (i == j) continue;
            IAtom jAtom = box.getLeafList().get(j);
            int jType = jAtom.getType().getIndex();
            Potential2Soft pij = ip[jType];
            if (pij == null) continue;
            if (bondingInfo.skipBondedPair(isPureAtoms, atom, jAtom)) continue;
            dr.Ev1Mv2(jAtom.getPosition(), atom.getPosition());
            boundary.nearestImage(dr);
            u += handleComputeOne(pij, zero, dr, zero, i, j);
        }
        return u;
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
            u += computeOneTruncationCorrection(atom.getLeafIndex());
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

    public void computeAllTruncationCorrection(double[] uCorrection, double[] duCorrection) {
        if (!doAllTruncationCorrection) {
            return;
        }
        double uCor = 0;
        double duCor = 0;
        for (int i = 0; i < atomCountByType.length; i++) {
            for (int j = i; j < atomCountByType.length; j++) {
                Potential2Soft p = pairPotentials[i][j];
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
