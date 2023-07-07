/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.nbr.cell;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.*;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.potential.BondingInfo;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.species.SpeciesManager;

import java.util.Arrays;

public class NeighborCellManager implements NeighborManager {
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
    protected final boolean rectangular;

    public NeighborCellManager(SpeciesManager sm, Box box, int cellRange, BondingInfo bondingInfo) {
        this.box = box;
        this.cellRange = cellRange;
        this.rectangular = box.getBoundary().isRectangular();
        boxHalf = box.getSpace().makeVector();
        numCells = new int[3];
        jump = new int[3];
        cellNextAtom = null;
        atomCell = cellOffsets = wrapMap = cellLastAtom = new int[0];
        this.isPureAtoms = sm.isPureAtoms();
        this.bondingInfo = bondingInfo;
        allCellOffsets = new int[0];
    }

    public Box getBox() {
        return box;
    }

    @Override
    public BondingInfo getBondingInfo() {
        return bondingInfo;
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
        return (az % numCells[2]) * numCells[0] * numCells[1]
                + (ay % numCells[1]) * numCells[0]
                + (ax % numCells[0]);
    }

    protected int wrappedIndex(int i, int nc) {
        int rv = i;
        while (rv < cellRange) rv += nc - 2 * cellRange;
        while (nc - rv <= cellRange) rv -= nc - 2 * cellRange;
        return rv;
    }

    @Override
    public NeighborIterator makeNeighborIterator() {
        return new NeighborIteratorCell(this, bondingInfo, isPureAtoms, box);
    }

    protected void addBoxListener() {
        box.getEventManager().addListener(new BoxEventListener() {
            @Override
            public void boxMoleculeAdded(BoxMoleculeEvent e) {
                if (cellNextAtom == null) return;
                int numAtoms = box.getLeafList().size();
                if (cellNextAtom.length < numAtoms) {
                    cellNextAtom = Arrays.copyOf(cellNextAtom, numAtoms);
                    atomCell = Arrays.copyOf(atomCell, numAtoms);
                }
                for (IAtom atom : e.getMolecule().getChildList()) {
                    int i = atom.getLeafIndex();
                    cellNextAtom[i] = -1;
                    atomCell[i] = -1;
                    updateAtom(atom);
                }
            }

            @Override
            public void boxMoleculeRemoved(BoxMoleculeEvent e) {
                for (IAtom atom : e.getMolecule().getChildList()) {
                    removeAtom(atom);
                }
            }

            @Override
            public void boxAtomLeafIndexChanged(BoxAtomIndexEvent e) {
                int oldIndex = e.getIndex();
                int newIndex = e.getAtom().getLeafIndex();
                moveAtomIndex(oldIndex, newIndex);
            }

            @Override
            public void boxNumberMolecules(BoxMoleculeCountEvent e) {
                if (cellNextAtom == null) return;
                int numAtoms = box.getLeafList().size() + e.getSpecies().getLeafAtomCount() * e.getCount();
                if (cellNextAtom.length < numAtoms) {
                    cellNextAtom = Arrays.copyOf(cellNextAtom, numAtoms);
                    atomCell = Arrays.copyOf(atomCell, numAtoms);
                }
            }
        });
    }

    public void init() {
        if (range == 0) {
            throw new RuntimeException("Need a range for init");
        }
        if (numCells[0] == 0) {
            addBoxListener();
        }
        int totalCells = 1;
        final Vector bs = box.getBoundary().getBoxSize();
        final boolean[] periodic = box.getBoundary().getPeriodicity();
        int D = periodic.length;
        // numCells is 3D even if we aren't, has 1 for extra components
        Arrays.fill(numCells, 1);
        Vector[] edgeVectors = null;
        int[] numBoxCells = null;
        double[] xmin = new double[D];
        double[] xmax = new double[D];
        if (rectangular) {
            double minCellSize = range / cellRange;
            for (int i = 0; i < D; i++) {
                // with cell lists, we can accommodate rc>L/2
                // we need the box to be at least the size of a cell.
                // when rc>L/2, we'll end up doing a lattice sum
                if (bs.getX(i) < minCellSize) {
                    throw new RuntimeException("box not big enough to accomodate even 1 cell (" + bs.getX(i) + " < " + minCellSize + ")");
                }
                // include cellRange of padding on each side of the box
                numCells[i] = ((int) Math.floor(bs.getX(i) / minCellSize));
                if (periodic[i]) numCells[i] += cellRange * 2;
                totalCells *= numCells[i];
            }
        } else {
            edgeVectors = new Vector[D];
            numBoxCells = new int[D];

            for (int i = 0; i < D; i++) {
                edgeVectors[i] = box.getBoundary().getEdgeVector(i);
            }

            // we assume that the edge vectors follow the pattern that the first (x)
            // is (Lxx,0,0) and the second (y) is (Lyx, Lyy, 0) and the last (z) is
            // (Lzx, Lyz, Lzz)

            // last dimension (z) first; lz = a Lz
            // min(lz) = minCellSize
            double minCellSize = range / cellRange;
            numBoxCells[D - 1] = numCells[D - 1] = (int) (edgeVectors[D - 1].getX(D - 1) / minCellSize);
            if (periodic[D - 1]) numCells[D - 1] += cellRange * 2;
            totalCells *= numCells[D - 1];

            xmin[D - 1] = 0;
            double lzz = edgeVectors[D - 1].getX(D - 1) / numBoxCells[D - 1];
            xmax[D - 1] = lzz;

            // now y
            // we need the minimum distance from ly=b Ly corner to the lz edge
            // the vector from the corner to the edge is ly - alpha lz and is perpendicular to lz
            // d^2 = (lyy - alpha lzy)^2 + alpha^2 lzz^2
            // d(d^2/dalpha) = 0
            // alpha = (lzy lyy) / (lzy^2 + lzz^2)
            // the vector is then
            // d = ly - (lzy lyy) / (lzy^2 + lzz^2) lz
            // it's squared magnitude is
            // d^2 = (lyy - (lzy lyy lzy) / (lzy^2 + lzz^2)) ^2 + ((lzy lyy lzz) / (lzy^2 + lzz^2))^2
            //     = lyy^2 ( (1 - lzy^2 / (lzy^2 + lzz^2))^2 + ((lzy lzz) / (lzy^2 + lzz^2))^2 )
            // we need this (d^2) to be at least minCellSize^2
            // lyy >= minCellSize / sqrt( (1 - lzy^2 / (lzy^2 + lzz^2))^2 + ((lzy lzz) / (lzy^2 + lzz^2))^2 )
            double lzy = edgeVectors[D - 1].getX(D - 2) / numBoxCells[D - 1];
            double facd = lzy * lzy + lzz * lzz;
            double fac1 = 1 - lzy * lzy / facd;
            double fac2 = lzy * lzz / facd;
            double minlyy = minCellSize / Math.sqrt(fac1 * fac1 + fac2 * fac2);
            numBoxCells[D - 2] = numCells[D - 2] = (int) (edgeVectors[D - 2].getX(D - 2) / minlyy);
            if (periodic[D - 2]) numCells[D - 2] += cellRange * 2;
            totalCells *= numCells[D - 2];

            double lyy = edgeVectors[D - 2].getX(D - 2) / numBoxCells[D - 2];
            xmin[D - 2] = Math.min(0, lzy);
            xmax[D - 2] = lyy;
            xmax[D - 2] += Math.max(0, lzy);

            if (D > 2) {
                // now x
                // we need the minimum distance from lx=c Lx corner to the ly, lz face
                // the vector from the corner to the face is lx - alpha lz - beta ly and is perpendicular to the face
                // (alpha can take a different value than it did for y above)
                // d^2 = (lxx - alpha lzx - beta lyx)^2 + (alpha lzy + beta lyy)^2 + alpha^2 lzz^2
                // d(d^2/dalpha) = 0       2 alpha lzx^2 - lxx lzx + beta lzx lyx + 2 alpha lzy^2 + 2 beta lzy lyy + 2 alpha lzz^2 = 0
                //                         (2 lzx^2 + 2 lzy^2 + 2 lzz^2) alpha + (lzx lyx + lzy lyy) beta = lxx lzx
                // d(d^2/dbeta) = 0        2 beta lyx^2 - lxx lyx + alpha lzx lyx + 2 beta lyy^2 + 2 alpha lzy lyy = 0
                //                         (lzx lyx + 2 lzy lyy) alpha + (2 lyx^2 + 2 lyy^2) beta = lxx lyx
                double lyx = edgeVectors[D - 2].getX(D - 3) / numBoxCells[D - 2];
                double lzx = edgeVectors[D - 1].getX(D - 3) / numBoxCells[D - 1];
                double a = 2 * lzx * lzx + 2 * lzy * lzy + 2 * lzz * lzz;
                double b = lzx * lyx + lzy * lyy;
                double c = lzx * lyx + 2 * lzy * lyy;
                double d = 2 * lyx * lyx + 2 * lyy * lyy;
                double det = a * d - b * c;
                double ai = d / det;
                double bi = -b / det;
                double ci = -c / det;
                double di = a / det;
                // alphaOlxx = alpha/lxx, betaOlxx = beta/lxx
                double alphaOlxx = ai * lzx + bi * lyx;
                double betaOlxx = ci * lzx + di * lyx;
                // d^2 = lxx^2 ( (1 - (alpha/lxx) lzx - (beta/lxx) lyx)^2 + ((alpha/lxx) lzy + (beta/lxx) lyy)^2 + (alpha/lxx)^2 lzz^2 ) >= minCellSize^2
                // lxx >= minCellSize / sqrt( (1 - (alpha/lxx) lzx - (beta/lxx) lyx)^2 + ((alpha/lxx) lzy + (beta/lxx) lyy)^2 + (alpha/lxx)^2 lzz^2 )
                fac1 = 1 - alphaOlxx * lzx - betaOlxx * lyx;
                fac2 = alphaOlxx * lzy + betaOlxx * lyy;
                double fac3 = alphaOlxx * lzz;
                double minlxx = minCellSize / Math.sqrt(fac1 * fac1 + fac2 * fac2 + fac3 * fac3);
                numBoxCells[D - 3] = numCells[D - 3] = (int) (edgeVectors[D - 3].getX(D - 3) / minlxx);
                if (periodic[D - 3]) numCells[D - 3] += cellRange * 2;

                xmin[D - 3] = 0;
                xmin[D - 3] += Math.min(0, lyx);
                xmin[D - 3] += Math.min(0, lzx);
                xmax[D - 3] = edgeVectors[D - 3].getX(D - 3) / numBoxCells[D - 3];
                xmax[D - 3] += Math.max(0, lyx);
                xmax[D - 3] += Math.max(0, lzx);
                totalCells *= numCells[D - 3];
            }
//            System.out.println("cells: "+numCells[0]+" "+numCells[1]+" "+numCells[2]);
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
            int izMax2 = D == 3 ? Math.min(icd2, cellRange * cellRange) : 0;
            for (int iz = 0; (iz2 = iz * iz) <= izMax2; iz++) {
                int izm1 = iz == 0 ? 0 : (Math.abs(iz) - 1);
                int izm1Sq = izm1 * izm1;
                int iy2Max = icd2 - iz2;
                int iyMax = D > 1 ? ((int) Math.sqrt(iy2Max + 0.001)) : 0;
                iyMax = Math.min(iyMax, cellRange);
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
                    ixMin = Math.max(ixMin, -cellRange);
                    ixMax = Math.min(ixMax, cellRange);
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
                        if (rectangular) {
                            if (ixm1Sq + iym1Sq + izm1Sq >= cellRange * cellRange) continue;
                        } else {
                            double dx = ix * edgeVectors[0].getX(0) / numBoxCells[0] + iy * edgeVectors[1].getX(0) / numBoxCells[1];
                            if (D > 2) dx += iz * edgeVectors[2].getX(0) / numBoxCells[2];
                            double min1 = xmin[0];
                            double max1 = xmax[0];
                            double min2 = min1 + dx;
                            double max2 = max1 + dx;
                            if ((max1 - min2) * (max2 - min1) > 0) dx = 0;
                            else if (dx > 0) dx = min2 - max1;
                            else dx = min1 - max2;
                            double dy = iy * edgeVectors[1].getX(1) / numBoxCells[1];
                            if (D > 2) dy += iz * edgeVectors[2].getX(1) / numBoxCells[2];
                            min1 = xmin[1];
                            max1 = xmax[1];
                            min2 = min1 + dy;
                            max2 = max1 + dy;
                            if ((max1 - min2) * (max2 - min1) > 0) dy = 0;
                            else if (dy > 0) dy = min2 - max1;
                            else dy = min1 - max2;
                            double dz = 0;
                            if (D > 2) {
                                dz = iz * edgeVectors[2].getX(2) / numBoxCells[2];
                                min1 = xmin[2];
                                max1 = xmax[2];
                                min2 = min1 + dz;
                                max2 = max1 + dz;
                                if ((max1 - min2) * (max2 - min1) > 0) dz = 0;
                                else if (dz > 0) dz = min2 - max1;
                                else dz = min1 - max2;
                            }
                            if (dx * dx + dy * dy + dz * dz > range * range) continue;
                        }
                        int mv = ix;
                        if (D > 0) {
                            if (D == 2) {
                                mv += iy * numCells[0];
                            } else if (D > 1) {
                                mv += (iy + iz * numCells[1]) * numCells[0];
                            }
                        }
//                        System.out.println(ix+" "+iy+" "+iz+" "+mv);

                        if (Math.abs(mv) > cellRange * numCells[0] * numCells[1] + cellRange * numCells[0] + cellRange) {
                            throw new RuntimeException("oops");
                        }
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
        int yboRange = D > 1 ? (periodic[1] ? (numCells[1] - cellRange - 1) / (numCells[1] - 2 * cellRange) : 0) : 0;
        int zboRange = D > 2 ? (periodic[2] ? (numCells[2] - cellRange - 1) / (numCells[2] - 2 * cellRange) : 0) : 0;
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
                    rawBoxOffsets[idx] = Vector.d(D);
                    if (rectangular) {
                        rawBoxOffsets[idx].setX(0, ix * bs.getX(0));
                        if (ny > 1) {
                            rawBoxOffsets[idx].setX(1, iy * bs.getX(1));
                        }
                        if (nz > 1) {
                            rawBoxOffsets[idx].setX(2, iz * bs.getX(2));
                        }
                    } else {
                        rawBoxOffsets[idx].Ea1Tv1(ix, edgeVectors[0]);
                        if (ny > 1) {
                            rawBoxOffsets[idx].PEa1Tv1(iy, edgeVectors[1]);
                        }
                        if (nz > 2) {
                            rawBoxOffsets[idx].PEa1Tv1(iz, edgeVectors[2]);
                        }
                    }
                }
            }
        }

        for (int ix = 0; ix < numCells[0]; ix++) {
            int x2 = periodic[0] ? wrappedIndex(ix, numCells[0]) : ix;
            int xbo = periodic[0] ? (ix - x2) / (numCells[0] - 2 * cellRange) : 0;
            for (int iy = 0; iy < numCells[1]; iy++) {
                int y2 = (D > 1 && periodic[1]) ? wrappedIndex(iy, numCells[1]) : iy;
                int ybo = (D > 1 && periodic[1]) ? (iy - y2) / (numCells[1] - 2 * cellRange) : 0;
                for (int iz = 0; iz < numCells[2]; iz++) {
                    int z2 = (D > 2 && periodic[2]) ? wrappedIndex(iz, numCells[2]) : iz;
                    int zbo = (D > 2 && periodic[2]) ? (iz - z2) / (numCells[2] - 2 * cellRange) : 0;
                    int iMap = iCell(ix, iy, iz);
                    int dCell = iCell(x2, y2, z2);
                    wrapMap[iMap] = dCell;
                    boxOffsets[iMap] = rawBoxOffsets[(xbo + xboRange) * ny * nz + (ybo + yboRange) * nz + (zbo + zboRange)];
                }
            }
        }

        if (allCellOffsets.length < 2 * numCellOffsets) {
            allCellOffsets = new int[2 * numCellOffsets];
        }
        System.arraycopy(cellOffsets, 0, allCellOffsets, 0, numCellOffsets);
        for (int i = 0; i < numCellOffsets; i++) {
            allCellOffsets[numCellOffsets + i] = -cellOffsets[i];
        }

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
        if (cellNextAtom == null || cellNextAtom.length < numAtoms) {
            cellNextAtom = new int[numAtoms];
            atomCell = new int[numAtoms];
        }

        boxHalf.Ea1Tv1(0.5, bs);
        jump[0] = 1;
        jump[1] = numCells[0];
        jump[2] = numCells[1] * numCells[0];
        for (int i = 0; i < cellLastAtom.length; i++) {
            cellLastAtom[i] = -1;
        }
        IAtomList atoms = box.getLeafList();
        Vector s = box.getSpace().makeVector();
        Tensor hInv = box.getBoundary().getHInv();
        for (int iAtom = 0; iAtom < numAtoms; iAtom++) {
            int cellNum = 0;
            IAtom atom = atoms.get(iAtom);
            Vector r = atom.getPosition();
            s.E(r);
            if (rectangular) {
                s.DE(bs);
            } else {
                hInv.transform(s);
            }
            for (int i = 0; i < r.getD(); i++) {
                double x = s.getX(i) + 0.5;
                int y = ((int) (cellRange + x * (numCells[i] - 2 * cellRange)));
                if (y == numCells[i] - cellRange) y--;
                else if (y == cellRange - 1) y++;
//                System.out.print(" "+y);
//                if (y < cellRange-1 || y > numCells[i] - cellRange) {
//                    throw new RuntimeException("oops");
//                }
                cellNum += y * jump[i];
            }
//            System.out.println();
            atomCell[iAtom] = cellNum;
             cellNextAtom[iAtom] = cellLastAtom[cellNum];
            cellLastAtom[cellNum] = iAtom;
        }
    }

    // only called from our box listener
    protected void removeAtom(IAtom atom) {
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
        Vector s = box.getSpace().makeVector();
        s.E(r);
        if (rectangular) {
            s.DE(bs);
        } else {
            Tensor hInv = box.getBoundary().getHInv();
            hInv.transform(s);
        }
        for (int i = 0; i < r.getD(); i++) {
            double x = s.getX(i) + 0.5;
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

    // only called from our box listener
    protected void moveAtomIndex(int oldIndex, int newIndex) {
        if (oldIndex == newIndex) return;
        int cell = atomCell[oldIndex];
        if (cellLastAtom[cell] == -1) {
            throw new RuntimeException("cell " + cell + " is empty, does not contain " + oldIndex);
        }
        // later we will try to remove this same "old" atom.  mark it as already removed.
        atomCell[oldIndex] = -1;
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
