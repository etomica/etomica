/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.nbr.cell;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.potential.BondingInfo;
import etomica.potential.Potential2Soft;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.simulation.Simulation;
import etomica.space.Vector;

import java.util.Arrays;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class NeighborCellManagerFasterer implements NeighborManager {
    protected final Box box;
    private final boolean isPureAtoms;
    private final BondingInfo bondingInfo;
    protected int cellRange;
    protected Vector[] boxOffsets;
    protected final Vector boxHalf;
    protected final int[] numCells, jump;
    protected double range;
    protected int[] wrapMap, cellLastAtom;
    protected int[] cellOffsets;
    protected int numCellOffsets;
    protected Vector[] rawBoxOffsets = new Vector[0];
    protected int[] cellNextAtom;
    protected int[] atomCell;
    public int[] allCellOffsets;

    public NeighborCellManagerFasterer(Simulation sim, Box box, int cellRange, BondingInfo bondingInfo) {
        this.box = box;
        this.cellRange = cellRange;
        boxHalf = box.getSpace().makeVector();
        numCells = new int[3];
        jump = new int[3];
        cellNextAtom = atomCell = cellOffsets = wrapMap = cellLastAtom = new int[0];
        this.isPureAtoms = sim.getSpeciesList().stream().allMatch(s -> s.getLeafAtomCount() == 1);
        this.bondingInfo = bondingInfo;
    }

    public int cellForCoord(Vector r) {
        int cellNum = 0;
        Vector bs = box.getBoundary().getBoxSize();
        boolean[] periodic = box.getBoundary().getPeriodicity();
        for (int i = 0; i < periodic.length; i++) {
            double x = (r.getX(i) + boxHalf.getX(i)) / bs.getX(i);
            if (periodic[i]) cellNum += ((int) (cellRange + x * (numCells[i] - 2 * cellRange))) * jump[i];
            else cellNum += ((int) (x * numCells[i])) * jump[i];
        }
        return cellNum;
    }

    public void setCellRange(int newCellRange) {
        cellRange = newCellRange;
    }

    public void setPotentialRange(double newRange) {
        range = newRange;
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

    private int iCell(int ax, int ay, int az) {
        return (az % numCells[2])
                + (ay % numCells[1]) * numCells[2]
                + (ax % numCells[0]) * numCells[1] * numCells[2];
    }

    protected int wrappedIndex(int i, int nc) {
        int rv = i;
        while (rv < cellRange) rv += nc - 2 * cellRange;
        while (nc - rv <= cellRange) rv -= nc - 2 * cellRange;
        return rv;
    }

    @Override
    public NeighborIterator makeNeighborIterator() {
        return new NeighborIteratorCellFasterer(this, bondingInfo, isPureAtoms, box);
    }

    public void init() {
        double minCellSize = range / cellRange;
        int totalCells = 1;
        final Vector bs = box.getBoundary().getBoxSize();
        final boolean[] periodic = box.getBoundary().getPeriodicity();
        for (int i = 0; i < 3; i++) {
            // with cell lists, we can accommodate rc>L/2
            // we need the box to be at least the size of a cell.
            // when rc>L/2, we'll end up doing a lattice sum
            if (bs.getX(i) < minCellSize) {
                System.err.println("box not big enough to accomodate even 1 cell (" + bs.getX(i) + " < " + minCellSize + ")");
                System.exit(1);
            }
            // include cellRange of padding on each side of the box
            numCells[i] = ((int) Math.floor(bs.getX(i) / minCellSize));
            if (periodic[i]) numCells[i] += cellRange * 2;
            totalCells *= numCells[i];
        }
        if (totalCells > cellLastAtom.length) {
            cellLastAtom = new int[totalCells];
            wrapMap = new int[totalCells];
            boxOffsets = new Vector[totalCells];
        }

        int lastCellCount = 0;
        numCellOffsets = 0;

        for (int icd = 1; icd <= cellRange + 5; icd++) {
            // we want all cells whose squared (index) distance is not more than icd
            // but exclude any cells handle in the previous passes
            int icd2 = icd * icd;
            int iz2;
            for (int iz = 0; (iz2 = iz * iz) <= icd2; iz++) {
                int izm1 = iz == 0 ? 0 : (Math.abs(iz) - 1);
                int izm1Sq = izm1 * izm1;
                int iy2Max = icd2 - iz2;
                int iyMax = (int) Math.sqrt(iy2Max + 0.001);
                for (int iy = -iyMax; iy <= iyMax; iy++) {
                    int iym1 = iy == 0 ? 0 : (Math.abs(iy) - 1);
                    int iym1Sq = iym1 * iym1;
                    int iy2 = iy * iy;
                    int ix2Max = iy2Max - iy2;
                    int ixMax = (int) Math.sqrt(ix2Max + 0.001);
                    int ix2Min = (icd - 1) * (icd - 1) - iz2 - iy2;
                    int ixMin;
                    if (ix2Min < 0) {
                        ix2Min = 0;
                        ixMin = -1;
                    } else {
                        ixMin = (int) Math.sqrt(ix2Min + 0.001);
                    }
                    for (int ix = -ixMax; ix <= ixMax; ix++) {
                        if (ix >= -ixMin && ix <= ixMin) {
                            ix = ixMin + 1;
                            if (ix > ixMax) break;
                        }
                        if (iz == 0) {
                            if (iy < 1) {
                                if (ix < 1) continue;
                            } else if (ix < 0) continue;
                        }
                        int ixm1 = ix == 0 ? 0 : (Math.abs(ix) - 1);
                        int ixm1Sq = ixm1 * ixm1;
                        if (ixm1Sq + iym1Sq + izm1Sq >= cellRange * cellRange) continue;
                        int mv = ix + (iy + iz * numCells[1]) * numCells[0];
                        // terrible? yes.
                        if (cellOffsets.length <= numCellOffsets)
                            cellOffsets = Arrays.copyOf(cellOffsets, numCellOffsets + 1);
                        cellOffsets[numCellOffsets] = mv;
                        numCellOffsets++;
                    }
                }
            }
            if (numCellOffsets == lastCellCount) break;
            lastCellCount = numCellOffsets;
        }

        int xboRange = periodic[0] ? (numCells[0] - cellRange - 1) / (numCells[0] - 2 * cellRange) : 0;
        int yboRange = periodic[1] ? (numCells[1] - cellRange - 1) / (numCells[1] - 2 * cellRange) : 0;
        int zboRange = periodic[2] ? (numCells[2] - cellRange - 1) / (numCells[2] - 2 * cellRange) : 0;
        int nx = (2 * xboRange + 1);
        int ny = (2 * yboRange + 1);
        int nz = (2 * zboRange + 1);
        if (rawBoxOffsets.length < nx * ny * nz) {
            rawBoxOffsets = new Vector[nx * ny * nz];
        }
        for (int ix = -xboRange; ix <= xboRange; ix++) {
            for (int iy = -yboRange; iy <= yboRange; iy++) {
                for (int iz = -zboRange; iz <= zboRange; iz++) {
                    int idx = (ix + xboRange) * ny * nz + (iy + yboRange) * nz + (iz + zboRange);
                    rawBoxOffsets[idx] = Vector.of(ix * bs.getX(0), iy * bs.getX(1), iz * bs.getX(2));
                }
            }
        }

        for (int ix = 0; ix < numCells[0]; ix++) {
            int x2 = periodic[0] ? wrappedIndex(ix, numCells[0]) : ix;
            int xbo = periodic[0] ? (ix - x2) / (numCells[0] - 2 * cellRange) : 0;
            for (int iy = 0; iy < numCells[1]; iy++) {
                int y2 = periodic[1] ? wrappedIndex(iy, numCells[1]) : iy;
                int ybo = periodic[1] ? (iy - y2) / (numCells[1] - 2 * cellRange) : 0;
                for (int iz = 0; iz < numCells[2]; iz++) {
                    int z2 = periodic[2] ? wrappedIndex(iz, numCells[2]) : iz;
                    int zbo = periodic[2] ? (iz - z2) / (numCells[2] - 2 * cellRange) : 0;
                    int iMap = iCell(ix, iy, iz);
                    int dCell = iCell(x2, y2, z2);
                    wrapMap[iMap] = dCell;
                    boxOffsets[iMap] = rawBoxOffsets[(xbo + xboRange) * ny * nz + (ybo + yboRange) * nz + (zbo + zboRange)];
                }
            }
        }

        this.allCellOffsets = Stream.of(
                IntStream.of(0),
                IntStream.of(cellOffsets),
                IntStream.of(cellOffsets).map(offset -> offset * -1)
        ).flatMapToInt(s -> s).toArray();

        assignCellAll();
    }

    public void assignCellAll() {
        Vector bs = box.getBoundary().getBoxSize();
        if (box.getBoundary().volume() == 0) {
            System.err.println("box has 0 volume, can't assign cells");
            return;
        }
        if (range == 0 || cellRange == 0) {
            System.err.println("range and cell range need to be non-zero");
            return;
        }

        final int numAtoms = box.getLeafList().size();
        if (cellNextAtom.length < numAtoms) {
            cellNextAtom = new int[numAtoms];
            atomCell = new int[numAtoms];
        }

        boxHalf.Ea1Tv1(0.5, bs);
        jump[0] = numCells[1] * numCells[2];
        jump[1] = numCells[2];
        jump[2] = 1;
        for (int i = 0; i < cellLastAtom.length; i++) {
            cellLastAtom[i] = -1;
        }
        IAtomList atoms = box.getLeafList();
        for (int iAtom = 0; iAtom < numAtoms; iAtom++) {
            int cellNum = 0;
            IAtom atom = atoms.get(iAtom);
            Vector r = atom.getPosition();
            for (int i = 0; i < 3; i++) {
                double x = (r.getX(i) + boxHalf.getX(i)) / bs.getX(i);
                int y = ((int) (cellRange + x * (numCells[i] - 2 * cellRange)));
                if (y == numCells[i] - cellRange) y--;
                else if (y == cellRange - 1) y++;
                cellNum += y * jump[i];
            }
            atomCell[iAtom] = cellNum;
            cellNextAtom[iAtom] = cellLastAtom[cellNum];
            cellLastAtom[cellNum] = iAtom;
        }
    }

    public void removeAtom(IAtom atom) {
        int iAtom = atom.getLeafIndex();
        int oldCell = atomCell[iAtom];
        if (oldCell > -1) {
            // delete from old cell
            int j = cellLastAtom[oldCell];
            if (j == iAtom) {
                cellLastAtom[oldCell] = cellNextAtom[iAtom];
            } else {
                while (cellNextAtom[j] != iAtom) {
                    j = cellNextAtom[j];
                }
                cellNextAtom[j] = cellNextAtom[iAtom];
            }
        }
    }

    public void updateAtom(IAtom atom) {
        int cellNum = 0;
        final Vector bs = box.getBoundary().getBoxSize();
        Vector r = atom.getPosition();
        for (int i = 0; i < 3; i++) {
            double x = (r.getX(i) + boxHalf.getX(i)) / bs.getX(i);
            int y = ((int) (cellRange + x * (numCells[i] - 2 * cellRange)));
            if (y == numCells[i] - cellRange) y--;
            else if (y == cellRange - 1) y++;
            cellNum += y * jump[i];
        }

        int iAtom = atom.getLeafIndex();
        int oldCell = atomCell[iAtom];
        // check if existing assignment is right
        if (cellNum == oldCell) return;
        if (oldCell > -1) {
            // delete from old cell
            int j = cellLastAtom[oldCell];
            if (j == iAtom) {
                cellLastAtom[oldCell] = cellNextAtom[iAtom];
            } else {
                while (cellNextAtom[j] != iAtom) {
                    j = cellNextAtom[j];
                }
                cellNextAtom[j] = cellNextAtom[iAtom];
            }
        }

        atomCell[iAtom] = cellNum;
        cellNextAtom[iAtom] = cellLastAtom[cellNum];
        cellLastAtom[cellNum] = iAtom;
    }

    public void moveAtomIndex(int oldIndex, int newIndex) {
        if (oldIndex == newIndex) return;
        int cell = atomCell[oldIndex];
        //printf("%d was in %d\n", oldIndex, cell);
        atomCell[newIndex] = cell;
        cellNextAtom[newIndex] = cellNextAtom[oldIndex];
        //printf("%d now points to %d\n", newIndex, cellNextAtom[newIndex]);
        if (cell == -1) return;
        int j = cellLastAtom[cell];
        if (j == oldIndex) {
            cellLastAtom[cell] = newIndex;
            return;
        }
        while (cellNextAtom[j] != oldIndex) j = cellNextAtom[j];
        cellNextAtom[j] = newIndex;
    }

    public Vector[] getBoxOffsets() {
        return boxOffsets;
    }

    public int[] getAtomCell() {
        return atomCell;
    }

    public int[] getCellNextAtom() {
        return cellNextAtom;
    }

    public int[] getCellOffsets() {
        return cellOffsets;
    }

    public int[] getWrapMap() {
        return wrapMap;
    }

    public int[] getCellLastAtom() {
        return cellLastAtom;
    }
}
