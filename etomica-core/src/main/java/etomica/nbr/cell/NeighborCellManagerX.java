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

public class NeighborCellManagerX extends NeighborCellManager {

    protected  int numCellsX;


    public NeighborCellManagerX(SpeciesManager sm, Box box, int cellRange, BondingInfo bondingInfo) {
        this (sm,box,cellRange,bondingInfo,false);
    }
    public NeighborCellManagerX(SpeciesManager sm, Box box, int cellRange, BondingInfo bondingInfo, boolean handleOutOfBox) {
        super(sm,box,cellRange,bondingInfo,handleOutOfBox);
    }










    private int iCell(int ax) {
        return
                 (ax );
    }








    public void init() {
        if (range == 0) {
            throw new RuntimeException("Need a range for init");
        }
        if (numCellsX == 0) {
            addBoxListener();
        }
        int totalCells = 1;
        final Vector bs = box.getBoundary().getBoxSize();
        final boolean[] periodic = box.getBoundary().getPeriodicity();
        int D = periodic.length;
        // numCells is 3D even if we aren't, has 1 for extra components
      numCellsX =1;
        Vector[] edgeVectors = null;
        int[] numBoxCells = null;
        double[] xmin = new double[D];
        double[] xmax = new double[D];
        if (rectangular) {
            double minCellSize = range / cellRange;
            for (int i = 0; i < 1; i++) {
                // with cell lists, we can accommodate rc>L/2
                // we need the box to be at least the size of a cell.
                // when rc>L/2, we'll end up doing a lattice sum
                if (bs.getX(i) < minCellSize) {
                    throw new RuntimeException("box not big enough to accomodate even 1 cell (" + bs.getX(i) + " < " + minCellSize + ")");
                }
                // include cellRange of padding on each side of the box
                numCellsX = ((int) Math.floor(bs.getX(i) / minCellSize));
                numCellsX += cellRange * 2;
                totalCells *= numCellsX;
            }
        }
//            System.out.println("cells: "+numCells[0]+" "+numCells[1]+" "+numCells[2]);

        if (handleOutOfBox){
            totalCells++;
        }
        if (totalCells > cellLastAtom.length) {
            cellLastAtom = new int[totalCells];
            wrapMap = new int[totalCells];
            boxOffsets = new Vector[totalCells];
        }

        int lastCellCount = 0;
        numCellOffsets = 0;


                    for (int ix = 1; ix <= cellRange; ix++) {


                        int mv = ix;

//                        System.out.println(ix+" "+iy+" "+iz+" "+mv);


                        // terrible? yes.
                        if (cellOffsets.length <= numCellOffsets)
                            cellOffsets = Arrays.copyOf(cellOffsets, numCellOffsets + 1);
                        cellOffsets[numCellOffsets] = mv;
                        numCellOffsets++;

            }




        int xboRange = periodic[0] ? (numCellsX - cellRange - 1) / (numCellsX - 2 * cellRange) : 0;

        int nx = (2 * xboRange + 1);

        if (rawBoxOffsets.length < nx ) {
            rawBoxOffsets = new Vector[nx ];
        }
        for (int ix = -xboRange; ix <= xboRange; ix++) {

                    int idx = (ix + xboRange) ;
                    rawBoxOffsets[idx] = Vector.d(D);
                    if (rectangular) {
                        rawBoxOffsets[idx].setX(0, ix * bs.getX(0));

                    } else {
                        rawBoxOffsets[idx].Ea1Tv1(ix, edgeVectors[0]);

                    }

        }

        for (int ix = 0; ix < numCellsX; ix++) {
            int x2 = periodic[0] ? wrappedIndex(ix, numCellsX) : ix;
            int xbo = periodic[0] ? (ix - x2) / (numCellsX - 2 * cellRange) : 0;

                    int iMap = iCell(ix);
                    int dCell = iCell(x2);
                    wrapMap[iMap] = dCell;
                    boxOffsets[iMap] = rawBoxOffsets[(xbo + xboRange) ];

        }
        if (handleOutOfBox){
            wrapMap[totalCells-1]=totalCells-1;
            boxOffsets[totalCells-1]= rawBoxOffsets[0];
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
            for (int i = 0; i < 1; i++) {
                double x = s.getX(i) + 0.5;
                int y = ((int) (cellRange + x * (numCellsX - 2 * cellRange)));
                if ( handleOutOfBox && (y<cellRange || y>= numCellsX -cellRange)){
                   cellNum=cellLastAtom.length-1;
                   break;
                }

                if (y == numCellsX - cellRange) y--;
                else if (y == cellRange - 1) y++;
//                System.out.print(" "+y);
//                if (y < cellRange-1 || y > numCells[i] - cellRange) {
//                    throw new RuntimeException("oops");
//                }
                cellNum += y ;

            }

//            System.out.println();
            atomCell[iAtom] = cellNum;
            cellNextAtom[iAtom] = cellLastAtom[cellNum];
            cellLastAtom[cellNum] = iAtom;
        }
    }

    // only called from our box listener


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
        for (int i = 0; i < 1; i++) {
            double x = s.getX(i) + 0.5;
            int y = cellRange +((int) ( x * (numCellsX - 2 * cellRange)));
            if ( handleOutOfBox && (y<cellRange || y>= numCellsX -cellRange)){
                cellNum=cellLastAtom.length-1;
                break;
            }
            if (y == numCellsX - cellRange) y--;
            else if (y == cellRange - 1) y++;
            cellNum += y ;

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

}
