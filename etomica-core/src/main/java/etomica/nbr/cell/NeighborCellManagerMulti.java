/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.nbr.cell;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.*;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.molecule.IMolecule;
import etomica.potential.BondingInfo;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.potential.compute.NeighborManagerCell;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.species.SpeciesManager;

import java.util.Arrays;

import static etomica.nbr.cell.NeighborIteratorCellMulti.doDebug;

public class NeighborCellManagerMulti implements NeighborManagerCell {
    protected final Box box;

    private final BondingInfo bondingInfo;
    protected int cellRange;
    public boolean doDebug= false;

    protected final Vector boxHalf;
    protected final int[] numCells, jump;
    protected double range;
    protected int[][] cellLastAtom;
    protected int[] cellOffsets;
    protected int numCellOffsets;

    protected int[] cellNextAtom;
    protected int[] atomCell;
    public int[] allCellOffsets;
    public int[] atomCellCoordinate;
    public int[][] moleculeCellOffsets;




    public NeighborCellManagerMulti(SpeciesManager sm, Box box, int cellRange, BondingInfo bondingInfo) {


        this.box = box;
        this.cellRange = cellRange;

        boxHalf = box.getSpace().makeVector();
        numCells = new int[3];
        jump = new int[3];
        cellNextAtom = null;
        atomCell = cellOffsets  = atomCellCoordinate =new int[0];
        cellLastAtom= new int[box.getMoleculeList().size()][0];
        moleculeCellOffsets= new int[box.getMoleculeList().size()][3];
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



    @Override
    public NeighborIterator makeNeighborIterator() {
       // return new NeighborIteratorCell(this, bondingInfo, isPureAtoms, box);
        return new NeighborIteratorCellMulti(this,bondingInfo,box);
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
                    atomCellCoordinate = Arrays.copyOf(atomCellCoordinate, numAtoms);
                }
                for (IAtom atom : e.getMolecule().getChildList()) {
                    int i = atom.getLeafIndex();
                    cellNextAtom[i] = -1;
                    atomCell[i] = -1;
                    atomCellCoordinate[i] = -1;
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
                throw new RuntimeException("Can't do that");
            }

            @Override
            public void boxNumberMolecules(BoxMoleculeCountEvent e) {
                if (cellNextAtom == null) return;
                int numAtoms = box.getLeafList().size() + e.getSpecies().getLeafAtomCount() * e.getCount();
                if (cellNextAtom.length < numAtoms) {
                    cellNextAtom = Arrays.copyOf(cellNextAtom, numAtoms);
                    atomCell = Arrays.copyOf(atomCell, numAtoms);
                    atomCellCoordinate = Arrays.copyOf(atomCellCoordinate, numAtoms);
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
            numCells[i] += cellRange * 2;
            totalCells *= numCells[i];
        }




        totalCells+=125;

        if (totalCells > cellLastAtom[0].length) {

            for(int i=0;i<cellLastAtom.length;i++) {
                cellLastAtom[i] = new int[totalCells];
            }
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

                        if (ixm1Sq + iym1Sq + izm1Sq >= cellRange * cellRange) continue;

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
            atomCellCoordinate = new int[numAtoms];
        }

        boxHalf.Ea1Tv1(0.5, bs);
        jump[0] = 1;
        jump[1] = numCells[0];
        jump[2] = numCells[1] * numCells[0];
        for (int i = 0; i < cellLastAtom.length; i++) {
            for(int j=0;j<cellLastAtom[i].length;j++){
                cellLastAtom[i][j] = -1;

            }

        }
        IAtomList atoms = box.getLeafList();
        Vector s = box.getSpace().makeVector();

        int yoffset=128- ((int) ((0.5 ) * (numCells[0] - 2 * cellRange)));

        for (int iAtom = 0; iAtom < numAtoms; iAtom++) {
            int cellNum = 0;
            IAtom atom = atoms.get(iAtom);
            Vector r = atom.getPosition();
            int m = atom.getParentGroup().getIndex();


            s.E(r);


            s.DE(bs);

            int z=0;
            int bigCellNum=0;
            int bigCellJump=1;
            for (int i = 0; i < r.getD(); i++) {
                double x = s.getX(i) + 0.5;

                int y = ((int) (cellRange + x * (numCells[i] - 2 * cellRange)));
                // ADD THIS DEBUG PRINT:
                if (doDebug && iAtom == 125) {
                    System.out.println("Atom 125: i=" + i + " x=" + x + " y=" + y);
                }
                if (atom.getIndex()==0) {
                    moleculeCellOffsets[m][i] = ((int) ((x-0.5 ) * (numCells[i] - 2 * cellRange)));
                    y -= moleculeCellOffsets[m][i];
                    cellNum += y * jump[i];

                    if  (false && moleculeCellOffsets[m][i] != 0 ) {
                        System.out.println(m + " " + i + " " + moleculeCellOffsets[m][i] + " " + y + " " + cellNum);
                    }
                }
                else {
                   y-=moleculeCellOffsets[m][i];
                    if (doDebug && iAtom == 1106) {
                        System.out.println("Atom 1106: i=" + i + " x=" + x + " y=" + y);
                    }
                    if (y < cellRange || y >= numCells[i] - cellRange) {
                        bigCellNum += (y < cellRange ? -1 : 1) * bigCellJump;


                    } else {

                        //if (y == numCells[i] - cellRange) y--;
                        //else if (y == cellRange - 1) y++;
//                System.out.print(" "+y);
//                if (y < cellRange-1 || y > numCells[i] - cellRange) {
//                    throw new RuntimeException("oops");
//                }
                        cellNum += y * jump[i];

                    }

                    bigCellJump *= 5;

                }
                z |= (y+yoffset) << (8 * i);
            }
            if (bigCellNum!=0){
                cellNum=cellLastAtom[m].length-63+bigCellNum;
               // z=-1;
            }
            if(cellNum<0){
              System.out.println(m+" "+ iAtom+" "+bigCellNum);
            }
//            System.out.println();
            atomCell[iAtom] = cellNum;
            atomCellCoordinate[iAtom]= z;
            cellNextAtom[iAtom] = cellLastAtom[m][cellNum];
            cellLastAtom[m][cellNum] = iAtom;
        }
    }

    // only called from our box listener
    protected void removeAtom(IAtom atom) {
        int iAtom = atom.getLeafIndex();
        int oldCell = atomCell[iAtom];
        if (oldCell > -1) {
            int m=atom.getParentGroup().getIndex();
            // delete from old cell
            int j = cellLastAtom[m][oldCell];
            if (j == iAtom) {
                cellLastAtom[m][oldCell] = cellNextAtom[iAtom];
            } else {
                while (cellNextAtom[j] != iAtom) {
                    j = cellNextAtom[j];
                }
                cellNextAtom[j] = cellNextAtom[iAtom];
            }
        }
    }

    public void updateAtom(IAtom atom) {
        if(atom.getIndex()==0){
            throw new RuntimeException("can not update core");
        }
        int cellNum = 0;
        final Vector bs = box.getBoundary().getBoxSize();
        Vector r = atom.getPosition();
        Vector s = box.getSpace().makeVector();
        IMolecule m = atom.getParentGroup();
        Vector r0 = m.getChildList().get(0).getPosition();
        s.Ev1Mv2(r,r0);
        int yoffset=128- ((int) ((0.5 ) * (numCells[0] - 2 * cellRange)));




        s.DE(bs);

        int z=0;
        int bigCellNum=0;
        int bigCellJump=1;

        for (int i = 0; i < r.getD(); i++) {
            double x = s.getX(i) + 0.5;
            int y = cellRange +((int) ( x * (numCells[i] - 2 * cellRange)));
            y-=moleculeCellOffsets[m.getIndex()][i];
            if(y<cellRange || y>= numCells[i]-cellRange){
                bigCellNum+=(y < cellRange?-1:1)*bigCellJump;

            }
            else {

                cellNum += y * jump[i];

            }
            z |= (y+yoffset) << (8 * i);

            bigCellJump*=5;

        }

        int mm=m.getIndex();
        if (bigCellNum!=0){
            cellNum=cellLastAtom[mm].length-63+bigCellNum;
            //z=-1;
        }

        int iAtom = atom.getLeafIndex();
        int oldCell = atomCell[iAtom];
        // check if existing assignment is right
        if (cellNum == oldCell) return;
        if (oldCell > -1) {
            // delete from old cell
            int j = cellLastAtom[mm][oldCell];
            if (j == iAtom) {
                cellLastAtom[mm][oldCell] = cellNextAtom[iAtom];
            } else {
                while (cellNextAtom[j] != iAtom) {
                    j = cellNextAtom[j];
                }
                cellNextAtom[j] = cellNextAtom[iAtom];
            }
        }

        atomCell[iAtom] = cellNum;
        atomCellCoordinate[iAtom]=z;
        cellNextAtom[iAtom] = cellLastAtom[mm][cellNum];
        cellLastAtom[mm][cellNum] = iAtom;
    }

    // only called from our box listener




    public int[] getAtomCell() {
        return atomCell;
    }

    public int[] getCellNextAtom() {
        return cellNextAtom;
    }

    public int[] getCellOffsets() {
        return cellOffsets;
    }



    public int getNumCellOffsets() {
        return numCellOffsets;
    }

    public int[][] getCellLastAtom() {
        return cellLastAtom;
    }


}
