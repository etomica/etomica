/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.nbr.cell;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.Potential2Soft;
import etomica.potential.PotentialMasterFasterer;
import etomica.simulation.Simulation;
import etomica.space.Vector;

public class PotentialMasterCellFasterer extends PotentialMasterFasterer {

    protected final NeighborCellManagerFasterer cellManager;
    protected int cellRange;

    public PotentialMasterCellFasterer(Simulation sim, Box box, int cellRange) {
        super(sim, box);
        cellManager = new NeighborCellManagerFasterer(box, cellRange);
        this.cellRange = cellRange;
    }

    public NeighborCellManagerFasterer getCellManager() {
        return cellManager;
    }

    @Override
    public void setPairPotential(AtomType atomType1, AtomType atomType2, Potential2Soft p12) {
        super.setPairPotential(atomType1, atomType2, p12);
        cellManager.setPotentialRange(getRange());
    }

    public double getRange() {
        double range = 0;
        for (int iType = 0; iType < pairPotentials.length; iType++) {
            for (int jType = iType; jType < pairPotentials.length; jType++) {
                if (pairPotentials[iType][jType] == null) continue;
                double rc = pairPotentials[iType][jType].getRange();
                if (rc > range) range = rc;
            }
        }
        return range;
    }

    public void init() {
        cellManager.init();
    }

    public void updateAtom(int iAtom) {
        cellManager.updateAtom(iAtom);
    }

    public double computeAll(boolean doForces) {
        zeroArrays(doForces);

        IAtomList atoms = box.getLeafList();
        double uTot = 0;
        Vector[] boxOffsets = cellManager.getBoxOffsets();
        int[] atomCell = cellManager.getAtomCell();
        int[] cellNextAtom = cellManager.getCellNextAtom();
        int[] cellOffsets = cellManager.getCellOffsets();
        int[] wrapMap = cellManager.getWrapMap();
        int[] cellLastAtom = cellManager.getCellLastAtom();
        for (int i = 0; i < atoms.size(); i++) {
            IAtom iAtom = atoms.get(i);
            Vector ri = positions.get(i);
            int iType = iAtom.getType().getIndex();
            Potential2Soft[] ip = pairPotentials[iType];
            int j = i;
            Vector jbo = boxOffsets[atomCell[i]];
            while ((j = cellNextAtom[j]) > -1) {
                IAtom jAtom = atoms.get(j);
                int jType = jAtom.getType().getIndex();
                Potential2Soft pij = ip[jType];
                if (pij == null) continue;
                if (skipBondedPair(isPureAtoms, iAtom, jAtom, bondedAtoms)) continue;
                uTot += handleComputeAll(doForces, i, j, ri, positions.get(j), jbo, pij);
            }
            int iCell = atomCell[i];
            for (int cellOffset : cellOffsets) {
                int jCell = iCell + cellOffset;
                jbo = boxOffsets[jCell];
                jCell = wrapMap[jCell];
                for (j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                    IAtom jAtom = atoms.get(j);
                    int jType = jAtom.getType().getIndex();
                    Potential2Soft pij = ip[jType];
                    if (pij == null) continue;
                    if (skipBondedPair(isPureAtoms, iAtom, jAtom, bondedAtoms)) continue;
                    uTot += handleComputeAll(doForces, i, j, ri, positions.get(j), jbo, pij);
                }
            }
        }
        if (!isPureAtoms) {
            uTot += computeAllBonds(doForces);
        }
        return uTot;
    }

    public double computeOne(IAtom iAtom) {
        int i = iAtom.getLeafIndex();
        Vector ri = positions.get(i);
        int iType = iAtom.getType().getIndex();
        IAtomList atoms = box.getLeafList();
        uAtomsChanged.clear();
//        uAtomsChanged.ensureCapacity(6000); // TODO
        uAtomsChanged.add(i);
        duAtom.clear();
//        duAtom.ensureCapacity(6000); // TODO
        duAtom.add(0);
        Potential2Soft[] ip = pairPotentials[iType];

        Vector[] boxOffsets = cellManager.getBoxOffsets();
        int[] atomCell = cellManager.getAtomCell();
        int[] cellNextAtom = cellManager.getCellNextAtom();
        int[] cellOffsets = cellManager.getCellOffsets();
        int[] wrapMap = cellManager.getWrapMap();
        int[] cellLastAtom = cellManager.getCellLastAtom();

        int iCell = atomCell[i];
        Vector jbo = boxOffsets[iCell];
        double u1 = 0;
        for (int j = cellLastAtom[iCell]; j > -1; j = cellNextAtom[j]) {
            if (j != i) {
                IAtom jAtom = atoms.get(j);
                int jType = jAtom.getType().getIndex();
                Potential2Soft pij = ip[jType];
                if (pij == null) continue;
                if (skipBondedPair(isPureAtoms, iAtom, jAtom, bondedAtoms)) continue;
                Vector rj = positions.get(j);
                u1 += handleComputeOne(pij, ri, rj, jbo, i, j);
            }
        }
        for (int cellOffset : cellOffsets) {
            int jCell = iCell + cellOffset;
            jbo = boxOffsets[jCell];
            jCell = wrapMap[jCell];
            for (int j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                IAtom jAtom = atoms.get(j);
                int jType = jAtom.getType().getIndex();
                Potential2Soft pij = ip[jType];
                if (pij == null) continue;
                if (skipBondedPair(isPureAtoms, iAtom, jAtom, bondedAtoms)) continue;
                u1 += handleComputeOne(pij, ri, positions.get(j), jbo, i, j);
            }

            // now down
            jCell = iCell - cellOffset;
            jbo = boxOffsets[jCell];
            jCell = wrapMap[jCell];
            for (int j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                IAtom jAtom = atoms.get(j);
                int jType = jAtom.getType().getIndex();
                Potential2Soft pij = ip[jType];
                if (pij == null) continue;
                if (skipBondedPair(isPureAtoms, iAtom, jAtom, bondedAtoms)) continue;
                u1 += handleComputeOne(pij, ri, positions.get(j), jbo, i, j);
            }
        }

        if (!isPureAtoms && !isOnlyRigidMolecules) {
            u1 += computeOneBonded(iAtom);
        }
//        System.out.println(uAtomsChanged.size());
        return u1;
    }
}
