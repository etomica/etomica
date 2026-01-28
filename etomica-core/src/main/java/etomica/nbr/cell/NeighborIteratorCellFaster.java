package etomica.nbr.cell;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.compute.NeighborIterator;
import etomica.space.Space;
import etomica.space.Vector;

public class NeighborIteratorCellFaster implements NeighborIterator {

    private final NeighborCellManager cellManager;


    private final Box box;
    private final Space space;
    private final boolean handleOutOfBox;



    public NeighborIteratorCellFaster(NeighborCellManager cellManager, Box box) {
        this.cellManager = cellManager;


        this.box = box;
        this.space = box.getSpace();
        this.handleOutOfBox = cellManager.isHandleOutOfBox();
    }



    @Override
    public void iterUpNeighbors(int iAtom, NeighborConsumer consumer) {

        IAtomList atoms = box.getLeafList();
        Vector[] boxOffsets = cellManager.getBoxOffsets();
        int[] atomCell = cellManager.getAtomCell();
        int[] cellNextAtom = cellManager.getCellNextAtom();
        int[] cellOffsets = cellManager.getCellOffsets();
        int[] wrapMap = cellManager.getWrapMap();
        int[] cellLastAtom = cellManager.getCellLastAtom();
        IAtom atom1 = atoms.get(iAtom);

        for (int j = cellNextAtom[iAtom]; j > -1; j = cellNextAtom[j]) {
            IAtom atom2 = atoms.get(j);

            Vector rij = space.makeVector();
            rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());

            consumer.accept(atom2, rij, 0);
        }

        int iCell = atomCell[iAtom];
        if (handleOutOfBox && iCell > wrapMap.length - 125) {

            for (int x = -1; x <= 1; x++) {
                for (int y = -1; y <= 1; y++) {
                    for (int z = 0; z <= 1; z++) {

                        int jCell = iCell + z * 25 + y * 5 + x;
                        if (jCell <= iCell) continue;
                        for (int j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                            IAtom atom2 = atoms.get(j);


                            Vector rij = space.makeVector();
                            rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());


                            consumer.accept(atom2, rij, 0);

                        }
                    }
                }
            }
            int c = iCell - wrapMap.length + 125;
            int x = c % 5;
            int y = ((c - x) / 5) % 5;
            int z = (((c - x) / 5) - y) % 5;
            int cr = cellManager.cellRange;
            int minX = cr, minY = cr, minZ = cr;
            int maxX = cellManager.numCells[0] - cr, maxY = cellManager.numCells[1] - cr, maxZ = cellManager.numCells[2] - cr;
            if (x == 1) maxX = minX + cr - 1;
            else if (x == 3) minX = maxX - cr + 1;
            if (y == 1) maxY = minY + cr - 1;
            else if (y == 3) minY = maxY - cr + 1;
            if (z == 1) maxZ = minZ + cr - 1;
            else if (z == 3) minZ = maxZ - cr + 1;
            for (int cX = minX; cX <= maxX; cX++) {
                for (int cY = minY; cY <= maxY; cY++) {
                    for (int cZ = minZ; cZ <= maxZ; cZ++) {
                        int jCell = cX + cY * cellManager.numCells[0] + cZ * cellManager.numCells[1] * cellManager.numCells[0];
                        for (int j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                            IAtom atom2 = atoms.get(j);


                            Vector rij = space.makeVector();
                            rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());


                            consumer.accept(atom2, rij, 0);

                        }
                    }
                }
            }
            return;


        }
        for (int ico = 0; ico < cellManager.numCellOffsets; ico++) {
            int cellOffset = cellOffsets[ico];
            int jCell = iCell + cellOffset;

            jCell = wrapMap[jCell];
            for (int j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                IAtom atom2 = atoms.get(j);


                Vector rij = space.makeVector();
                rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());


                consumer.accept(atom2, rij, 0);

            }
        }
    }

    @Override
    public void iterDownNeighbors(int iAtom, NeighborConsumer consumer) {

        IAtomList atoms = box.getLeafList();
        Vector[] boxOffsets = cellManager.getBoxOffsets();
        int[] atomCell = cellManager.getAtomCell();
        int[] cellNextAtom = cellManager.getCellNextAtom();
        int[] cellOffsets = cellManager.getCellOffsets();
        int[] wrapMap = cellManager.getWrapMap();
        int[] cellLastAtom = cellManager.getCellLastAtom();
        IAtom atom1 = atoms.get(iAtom);

        int iCell = atomCell[iAtom];
        for (int j = cellLastAtom[iCell]; j != iAtom; j = cellNextAtom[j]) {
            IAtom atom2 = atoms.get(j);

            Vector rij = space.makeVector();
            rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());

            consumer.accept(atom2, rij, 0);
        }
        if (handleOutOfBox && iCell==wrapMap.length-1){
            for (int jCell = 0; jCell < wrapMap.length-1; jCell++) {



                for (int j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                    IAtom atom2 = atoms.get(j);


                    Vector rij = space.makeVector();
                    rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());

                    consumer.accept(atom2, rij, 0);
                }
            }
            return;
        }

        for (int ico = 0; ico < cellManager.numCellOffsets; ico++) {
            int cellOffset = cellOffsets[ico];
            int jCell = iCell - cellOffset;
            Vector jbo = boxOffsets[jCell];
            jCell = wrapMap[jCell];
            for (int j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                IAtom atom2 = atoms.get(j);


                Vector rij = space.makeVector();
                rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                rij.PE(jbo);
                consumer.accept(atom2, rij, 0);
            }
        }
    }

    @Override
    public double iterAndSumAllNeighbors(IAtom atom1, SuperNbrConsumer consumer) {

        IAtomList atoms = box.getLeafList();
        Vector[] boxOffsets = cellManager.getBoxOffsets();
        int[] atomCell = cellManager.getAtomCell();
        int[] cellNextAtom = cellManager.getCellNextAtom();
//        int[] cellOffsets = cellManager.getCellOffsets();
        int[] allCellOffsets = cellManager.allCellOffsets;
        int[] wrapMap = cellManager.getWrapMap();
        int[] cellLastAtom = cellManager.getCellLastAtom();
        int iAtom = atom1.getLeafIndex();
        int iCell = atomCell[iAtom];
        double sum = 0;

//        for (int cellOffset : allCellOffsets) {
//            int jCell = iCell + cellOffset;
//            Vector jbo = boxOffsets[jCell];
//            jCell = wrapMap[jCell];
//            for (int j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
//                if (j == iAtom) { continue; }
//                IAtom atom2 = atoms.get(j);
//                if (bondingInfo.skipBondedPair(isPureAtoms, atom1, atom2)) {
//                    continue;
//                }
//
//                Vector rij = space.makeVector();
//                rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
//                rij.PE(jbo);
//                sum += consumer.accept(atom1, atom2, rij);
//            }
//        }
//        return sum;


        for (int j = cellLastAtom[iCell]; j > -1; j = cellNextAtom[j]) {
            if (j == iAtom) {
                continue;
            }
            IAtom atom2 = atoms.get(j);

            Vector rij = space.makeVector();
            rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());

            sum += consumer.accept(atom1, atom2, rij, 0);
        }

        for (int ico = 0; ico < 2 * cellManager.numCellOffsets; ico++) {
            int cellOffset = allCellOffsets[ico];
            int jCell = iCell + cellOffset;
            Vector jbo = boxOffsets[jCell];
            jCell = wrapMap[jCell];
            for (int j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                IAtom atom2 = atoms.get(j);


                Vector rij = space.makeVector();
                rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                rij.PE(jbo);

                sum += consumer.accept(atom1, atom2, rij, 0);
            }

//            jCell = iCell - cellOffset;
//            jbo = boxOffsets[jCell];
//            jCell = wrapMap[jCell];
//            for (int j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
//                IAtom atom2 = atoms.get(j);
//                if (bondingInfo.skipBondedPair(isPureAtoms, atom1, atom2)) {
//                    continue;
//                }
//
//                Vector rij = space.makeVector();
//                rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
//                rij.PE(jbo);
//                sum += consumer.accept(atom1, atom2, rij);
//            }
        }
        return sum;
    }

    @Override
    public void iterAllNeighbors(int iAtom, NeighborConsumer consumer) {

        IAtomList atoms = box.getLeafList();
        Vector[] boxOffsets = cellManager.getBoxOffsets();
        int[] atomCell = cellManager.getAtomCell();
        int[] cellNextAtom = cellManager.getCellNextAtom();
        int[] allCellOffsets = cellManager.allCellOffsets;
        int[] wrapMap = cellManager.getWrapMap();
        int[] cellLastAtom = cellManager.getCellLastAtom();
        IAtom atom1 = atoms.get(iAtom);
        int iCell = atomCell[iAtom];

        for (int j = cellLastAtom[iCell]; j > -1; j = cellNextAtom[j]) {
            if (j == iAtom) {
                continue;
            }
            IAtom atom2 = atoms.get(j);

            Vector rij = space.makeVector();
            rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());

            consumer.accept(atom2, rij, 0);
        }
        if (handleOutOfBox && iCell==wrapMap.length-1){
            for (int jCell = 0; jCell < wrapMap.length-1; jCell++) {



                for (int j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                    IAtom atom2 = atoms.get(j);


                    Vector rij = space.makeVector();
                    rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());

                    consumer.accept(atom2, rij, 0);
                }
            }
            return;
        }

        for (int ico = 0; ico < 2 * cellManager.numCellOffsets; ico++) {
            int cellOffset = allCellOffsets[ico];
            int jCell = iCell + cellOffset;
            Vector jbo = boxOffsets[jCell];
            jCell = wrapMap[jCell];
            for (int j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                IAtom atom2 = atoms.get(j);


                Vector rij = space.makeVector();
                rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                rij.PE(jbo);

                consumer.accept(atom2, rij, 0);
            }
        }
        if (handleOutOfBox){
            int jCell=wrapMap.length-1;
            for (int j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                IAtom atom2 = atoms.get(j);


                Vector rij = space.makeVector();
                rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());

                consumer.accept(atom2, rij, 0);

            }
        }
    }

}
