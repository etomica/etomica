/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.BondingInfo;
import etomica.potential.IPotential2;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.species.SpeciesManager;

/**
 * Subclass that handles path integrals by skipping any atom pair with a different child index.
 */
public class PotentialMasterCellPI extends PotentialMasterCell {

    public PotentialMasterCellPI(SpeciesManager sm, Box box, int cellRange, BondingInfo bondingInfo) {
        super(sm, box, cellRange, bondingInfo);
    }

    public double computeAll(boolean doForces, PotentialCallback pc) {
        double[] uAtomOld = new double[uAtom.length];
        boolean debug = false;
        if (debug) System.arraycopy(uAtom, 0, uAtomOld, 0, uAtom.length);
        zeroArrays(doForces);

        IAtomList atoms = box.getLeafList();
        double uTot = 0;
        Vector[] boxOffsets = cellManager.getBoxOffsets();
        int[] atomCell = cellManager.getAtomCell();
        int[] cellNextAtom = cellManager.getCellNextAtom();
        int[] cellOffsets = cellManager.getCellOffsets();
        int[] wrapMap = cellManager.getWrapMap();
        int[] cellLastAtom = cellManager.getCellLastAtom();
        int numCellOffsets = cellManager.getNumCellOffsets();
        long t1 = System.nanoTime();
        for (int i = 0; i < atoms.size(); i++) {
            IAtom iAtom = atoms.get(i);
            Vector ri = iAtom.getPosition();
            int iType = iAtom.getType().getIndex();
            IPotential2[] ip = pairPotentials[iType];
            int j = i;
            Vector jbo = boxOffsets[atomCell[i]];
            while ((j = cellNextAtom[j]) > -1) {
                IAtom jAtom = atoms.get(j);
                if (iAtom.getIndex() != jAtom.getIndex()) continue;
                int jType = jAtom.getType().getIndex();
                IPotential2 pij = ip[jType];
                if (pij == null) continue;
                boolean skipIntra = bondingInfo.skipBondedPair(isPureAtoms, iAtom, jAtom);
                uTot += handleComputeAll(doForces, i, j, ri, jAtom.getPosition(), jbo, pij, pc, skipIntra);
            }
            int iCell = atomCell[i];
            for (int k=0; k<numCellOffsets; k++) {
                int cellOffset = cellOffsets[k];
                int jCell = iCell + cellOffset;
                jbo = boxOffsets[jCell];
                jCell = wrapMap[jCell];
                for (j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                    IAtom jAtom = atoms.get(j);
                    if (iAtom.getIndex() != jAtom.getIndex()) continue;
                    int jType = jAtom.getType().getIndex();
                    IPotential2 pij = ip[jType];
                    if (pij == null) continue;
                    boolean skipIntra = bondingInfo.skipBondedPair(isPureAtoms, iAtom, jAtom);
                    uTot += handleComputeAll(doForces, i, j, ri, jAtom.getPosition(), jbo, pij, pc, skipIntra);
                }
            }
        }
        tAll += System.nanoTime() - t1;
        double[] uCorrection = new double[1];
        double[] duCorrection = new double[1];
        this.computeAllTruncationCorrection(uCorrection, duCorrection);
        uTot += uCorrection[0];
        virialTot += duCorrection[0];
        if (doForces && !isPureAtoms) {
            virialTot += PotentialCompute.computeVirialIntramolecular(forces, box);
        }

        if (debug && uAtom.length == uAtomOld.length) {
            boolean success = true;
            for (int i = 0; i < uAtom.length; i++) {
                if (Math.abs(uAtom[i] - uAtomOld[i]) > 1e-9) {
                    System.out.println("uAtom diff " + i + " " + uAtom[i] + " " + uAtomOld[i]);
                    success = false;
                }
            }
            if (!success) throw new RuntimeException("oops");
            System.out.println("success!");
        }
        energyTot = uTot;
        return uTot;
    }

    protected double computeOneInternal(IAtom iAtom, int startExcludeIdx, IAtom... excludedAtoms) {
        int iType = iAtom.getType().getIndex();
        int i = iAtom.getLeafIndex();
        Vector ri = iAtom.getPosition();
        IAtomList atoms = box.getLeafList();
        IPotential2[] ip = pairPotentials[iType];

        Vector[] boxOffsets = cellManager.getBoxOffsets();
        int[] atomCell = cellManager.getAtomCell();
        int[] cellNextAtom = cellManager.getCellNextAtom();
        int[] cellOffsets = cellManager.getCellOffsets();
        int[] wrapMap = cellManager.getWrapMap();
        int[] cellLastAtom = cellManager.getCellLastAtom();

        int iCell = atomCell[i];
        Vector jbo = boxOffsets[iCell];
        double u1 = 0;
        long t1 = System.nanoTime();
        for (int j = cellLastAtom[iCell]; j > -1; j = cellNextAtom[j]) {
            if (j != i) {
                IAtom jAtom = atoms.get(j);
                if (iAtom.getIndex() != jAtom.getIndex()) continue;
                int jType = jAtom.getType().getIndex();
                IPotential2 pij = ip[jType];
                if (pij == null) continue;
                if (arrayContains(jAtom, startExcludeIdx, excludedAtoms)) continue;
                boolean skipIntra = bondingInfo.skipBondedPair(isPureAtoms, iAtom, jAtom);
                Vector rj = jAtom.getPosition();
                u1 += handleComputeOne(pij, ri, rj, jbo, i, j, skipIntra);
            }
        }
        for (int cellOffset : cellOffsets) {
            int jCell = iCell + cellOffset;
            jbo = boxOffsets[jCell];
            jCell = wrapMap[jCell];
            for (int j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                IAtom jAtom = atoms.get(j);
                if (iAtom.getIndex() != jAtom.getIndex()) continue;
                int jType = jAtom.getType().getIndex();
                IPotential2 pij = ip[jType];
                if (pij == null) continue;
                if (arrayContains(jAtom, startExcludeIdx, excludedAtoms)) continue;
                boolean skipIntra = bondingInfo.skipBondedPair(isPureAtoms, iAtom, jAtom);
                u1 += handleComputeOne(pij, ri, jAtom.getPosition(), jbo, i, j, skipIntra);
            }

            // now down
            jCell = iCell - cellOffset;
            jbo = boxOffsets[jCell];
            jCell = wrapMap[jCell];
            for (int j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                IAtom jAtom = atoms.get(j);
                if (iAtom.getIndex() != jAtom.getIndex()) continue;
                int jType = jAtom.getType().getIndex();
                IPotential2 pij = ip[jType];
                if (pij == null) continue;
                if (arrayContains(jAtom, startExcludeIdx, excludedAtoms)) continue;
                boolean skipIntra = bondingInfo.skipBondedPair(isPureAtoms, iAtom, jAtom);
                u1 += handleComputeOne(pij, ri, jAtom.getPosition(), jbo, i, j, skipIntra);
            }
        }
        tMC += System.nanoTime() - t1;
        return u1;
    }
}
