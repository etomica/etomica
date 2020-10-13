/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.box.storage.IntStorage;
import etomica.box.storage.Tokens;
import etomica.box.storage.VectorStorage;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;

import java.util.*;

public class PotentialMasterFasterer {
    private final List<ISpecies> speciesList;
    protected final Potential2Soft[][] pairPotentials;
    protected final Box box;
    protected final Vector dr;
    protected double[] uAtom;
    protected final DoubleArrayList duAtom;
    protected final IntArrayList uAtomsChanged;
    protected Vector zero;
    protected double virialTot;
//    protected Vector[] forces;
    protected final VectorStorage forces;
    protected final Space space;

    protected final boolean isPureAtoms;
    protected boolean isOnlyRigidMolecules = true;
    protected int[][][] bondedAtoms;
    protected Map<Potential2Soft, List<int[]>>[] bondedPairs;
    protected Map<Potential2Soft, int[][]>[] bondedPartners;
    protected final VectorStorage positions;
    protected final IntStorage atomTypeIds;

    public PotentialMasterFasterer(Simulation sim, Box box) {
        space = box.getSpace();
        this.speciesList = sim.getSpeciesList();
        ISpecies species = sim.getSpecies(sim.getSpeciesCount() - 1);
        int lastTypeIndex = species.getAtomType(species.getUniqueAtomTypeCount() - 1).getIndex();
        pairPotentials = new Potential2Soft[lastTypeIndex + 1][lastTypeIndex + 1];
        this.box = box;
        positions = box.getAtomStorage(Tokens.POSITION);
        atomTypeIds = box.getAtomStorage(Tokens.ATOM_TYPE_ID);
        dr = sim.getSpace().makeVector();

        uAtom = new double[box.getLeafList().size()];
        uAtomsChanged = new IntArrayList(16);
        duAtom = new DoubleArrayList(16);
        zero = box.getSpace().makeVector();
        forces = box.getAtomStorage(Tokens.vectorsDefault());

        isPureAtoms = speciesList.stream().allMatch(s -> s.getLeafAtomCount() == 1);
        bondedAtoms = new int[sim.getSpeciesCount()][][];
        bondedPairs = new HashMap[sim.getSpeciesCount()];
        bondedPartners = new HashMap[sim.getSpeciesCount()];
        for (int i = 0; i < bondedPairs.length; i++) {
            bondedPairs[i] = new HashMap<>();
            bondedPartners[i] = new HashMap<>();
        }

    }

    public void init() {
    }

    public VectorStorage getForces() {
        return forces;
    }

    public double getLastVirial() {
        return virialTot;
    }

    public double getOldEnergy() {
        double uTot = 0;
        for (double iuAtom : uAtom) {
            uTot += iuAtom;
        }
        return uTot;
    }

    public void setPairPotential(AtomType atomType1, AtomType atomType2, Potential2Soft p12) {
        pairPotentials[atomType1.getIndex()][atomType2.getIndex()] = p12;
    }

    public void setBondingPotential(ISpecies species, Potential2Soft potential, List<int[]> bondedIndices) {
        isOnlyRigidMolecules = false;
        int speciesIndex = species.getIndex();
        if (bondedAtoms[speciesIndex] == null) {
            bondedAtoms[speciesIndex] = new int[species.getLeafAtomCount()][0];
        }

        int[][] partners = new int[species.getLeafAtomCount()][0];

        int[][] speciesAtoms = bondedAtoms[speciesIndex];
        for (int[] indices : bondedIndices) {
            int iAtom = Math.min(indices[0], indices[1]);
            int jAtom = Math.max(indices[0], indices[1]);
            speciesAtoms[iAtom] = etomica.util.Arrays.addInt(speciesAtoms[iAtom], jAtom);
            partners[iAtom] = etomica.util.Arrays.addInt(partners[iAtom], jAtom);
            partners[jAtom] = etomica.util.Arrays.addInt(partners[jAtom], iAtom);
        }
        if (bondedPairs[speciesIndex].containsKey(potential)) {
            throw new RuntimeException("Attempting to add the same bonding potential twice");
        }
        bondedPairs[speciesIndex].put(potential, new ArrayList<>(bondedIndices));
        bondedPartners[speciesIndex].put(potential, partners);

    }

    public void updateAtom(int iAtom) {

    }

    protected double handleComputeAll(boolean doForces, int iAtom, int jAtom, Vector ri, Vector rj, Vector jbo, Potential2Soft pij) {
        Vector dr = space.makeVector();
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
            forces.get(iAtom).PE(dr);
            forces.get(iAtom).ME(dr);
        }
        return uij;
    }

    public double computeAllBonds(boolean doForces) {
        double[] uTot = {0};
        for (int i = 0; i < speciesList.size(); i++) {
            Map<Potential2Soft, List<int[]>> potentials = bondedPairs[i];
            IMoleculeList molecules = box.getMoleculeList(speciesList.get(i));
            potentials.forEach((potential, pairs) -> {
                for (IMolecule molecule : molecules) {
                    for (int[] pair : pairs) {
                        IAtom iAtom = molecule.getChildList().get(pair[0]);
                        IAtom jAtom = molecule.getChildList().get(pair[1]);
                        uTot[0] += handleOneBondPair(doForces, iAtom.getLeafIndex(), jAtom.getLeafIndex(), potential);
                    }
                }
            });
        }
        return uTot[0];
    }

    protected double computeOneBonded(IAtom atom) {
        double[] uTot = {0};
        IMolecule parentMolecule = atom.getParentGroup();
        bondedPartners[atom.getParentGroup().getType().getIndex()].forEach((potential, partners) -> {
            for (int partnerIdx : partners[atom.getIndex()]) {
                IAtom jAtom = parentMolecule.getChildList().get(partnerIdx);
                uTot[0] += handleComputeOneBonded(potential, atom.getLeafIndex(), jAtom.getLeafIndex());
            }
        });
        return uTot[0];
    }

    protected double handleComputeOneBonded(Potential2Soft pij, int iAtom, int jAtom) {
        Vector ri = positions.get(iAtom);
        Vector rj = positions.get(jAtom);
        dr.Ev1Mv2(rj, ri);
        box.getBoundary().nearestImage(dr);
        double uij = pij.u(dr.squared());
        if (uij == 0) return 0;
        duAtom.plusEquals(0, 0.5 * uij);
        duAtom.add(0.5 * uij);
        uAtomsChanged.add(jAtom);
        return uij;
    }

    private double handleOneBondPair(boolean doForces, int iAtom, int jAtom, Potential2Soft potential) {
        Vector ri = positions.get(iAtom);
        Vector rj = positions.get(jAtom);
        dr.Ev1Mv2(rj, ri);
        box.getBoundary().nearestImage(dr);
        double r2 = dr.squared();
        double[] u = {0};
        double[] du = {0};
        potential.udu(r2, u, du);
        double uij = u[0];
        if (uij == 0) return 0;
        uAtom[iAtom] += 0.5 * uij;
        uAtom[jAtom] += 0.5 * uij;

        if (doForces) {
            double duij = du[0];
            dr.TE(duij / r2);
            forces.get(iAtom).PE(dr);
            forces.get(jAtom).ME(dr);
        }
        return uij;
    }

    protected static boolean skipBondedPair(boolean isPureAtoms, IAtom iAtom, IAtom jAtom, int[][][] bondedAtoms) {
        if (!isPureAtoms && iAtom.getParentGroup() == jAtom.getParentGroup()) {
            // ensure i < j
            if (iAtom.getLeafIndex() > jAtom.getLeafIndex()) {
                IAtom tmp = iAtom;
                iAtom = jAtom;
                jAtom = tmp;
            }
            int species = iAtom.getParentGroup().getType().getIndex();
            int iChildIndex = iAtom.getIndex();
            int jChildIndex = jAtom.getIndex();
            int[] iBondedAtoms = bondedAtoms[species][iChildIndex];
            for (int iBondedAtom : iBondedAtoms) {
                if (iBondedAtom == jChildIndex) return true;
            }
        }
        return false;
    }

    protected final void zeroArrays(boolean doForces) {
        virialTot = 0;

        int numAtoms = box.getLeafList().size();
        if (numAtoms > uAtom.length) {
            uAtom = new double[numAtoms];
        }
        for (int i = 0; i < numAtoms; i++) {
            uAtom[i] = 0;
            if (doForces) forces.get(i).E(0);
        }
    }

    public double computeAll(boolean doForces) {
        zeroArrays(doForces);

        IAtomList atoms = box.getLeafList();
        double u = 0;
        Boundary boundary = box.getBoundary();
        for (int i = 0; i < atoms.size(); i++) {
            IAtom iAtom = atoms.get(i);
            Vector ri = positions.get(i);
            int iType = iAtom.getType().getIndex();
            Potential2Soft[] ip = pairPotentials[iType];
            for (int j = i + 1; j < atoms.size(); j++) {
                IAtom jAtom = atoms.get(j);
                int jType = jAtom.getType().getIndex();
                Potential2Soft pij = ip[jType];
                if (pij == null) continue;

                if (skipBondedPair(isPureAtoms, iAtom, jAtom, bondedAtoms)) continue;

                dr.Ev1Mv2(positions.get(j), ri);
                boundary.nearestImage(dr);
                u += handleComputeAll(doForces, i, j, zero, dr, zero, pij);
            }
        }

        if (!isPureAtoms) {
            u += computeAllBonds(doForces);
        }

        return u;
    }

    public double computeOneOld(int iAtom) {
        return uAtom[iAtom] * 2;
    }

    protected double handleComputeOne(Potential2Soft pij, Vector ri, Vector rj, Vector jbo, int iAtom, int jAtom) {
        dr.Ev1Mv2(rj, ri);
        dr.PE(jbo);
        double uij = pij.u(dr.squared());
        if (uij == 0) return 0;

        duAtom.plusEquals(0, 0.5 * uij);
        duAtom.add(0.5 * uij);
        uAtomsChanged.add(jAtom);
        return uij;
    }

    public double computeOne(int iAtom) {
        IAtom a = box.getLeafList().get(iAtom);
        Vector ri = positions.get(iAtom);
        int iType = a.getType().getIndex();
        IAtomList atoms = box.getLeafList();
        uAtomsChanged.clear();
        uAtomsChanged.ensureCapacity(atoms.size());
        uAtomsChanged.add(iAtom);
        duAtom.clear();
        duAtom.ensureCapacity(atoms.size());
        duAtom.add(0);

        Potential2Soft[] ip = pairPotentials[iType];
        double u = 0;
        Boundary boundary = box.getBoundary();
        int s = atoms.size();
        for (int j = 0; j < s; j++) {
            if (iAtom == j) continue;
            IAtom jAtom = atoms.get(j);
            int jType = jAtom.getType().getIndex();
            Potential2Soft pij = ip[jType];
            if (pij == null) continue;
            if (skipBondedPair(isPureAtoms, a, jAtom, bondedAtoms)) continue;
            dr.Ev1Mv2(positions.get(j), ri);
            boundary.nearestImage(dr);
            u += handleComputeOne(pij, zero, dr, zero, iAtom, j);
        }
        if (!isPureAtoms && !isOnlyRigidMolecules) {
            u += computeOneBonded(a);
        }
        return u;
    }

    public void processAtomU(double fac) {
        for (int j = 0; j < uAtomsChanged.size(); j++) {
            int jj = uAtomsChanged.getInt(j);
            uAtom[jj] += fac * duAtom.getDouble(j);
        }
    }
}
