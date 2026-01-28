package etomica.nbr.cell;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.compute.NeighborIterator;
import etomica.space.Space;
import etomica.space.Vector;

public class NeighborIteratorCellFaster1 implements NeighborIterator {

    private final NeighborCellManager1 cellManager;


    private final Box box;
    private final Space space;
    private final boolean handleOutOfBox;



    public NeighborIteratorCellFaster1(NeighborCellManager1 cellManager, Box box) {
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
        if (handleOutOfBox && iCell==wrapMap.length-1){


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
        if (handleOutOfBox ){
            int jCell=wrapMap.length-1;
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
