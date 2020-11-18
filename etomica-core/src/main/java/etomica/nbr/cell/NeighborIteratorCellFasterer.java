package etomica.nbr.cell;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.BondingInfo;
import etomica.potential.compute.NeighborIterator;
import etomica.space.Space;
import etomica.space.Vector;

public class NeighborIteratorCellFasterer implements NeighborIterator {

    private final NeighborCellManagerFasterer cellManager;
    private final BondingInfo bondingInfo;
    private final boolean isPureAtoms;
    private final Box box;
    private final Space space;

    public NeighborIteratorCellFasterer(NeighborCellManagerFasterer cellManager, BondingInfo bondingInfo, boolean isPureAtoms, Box box) {
        this.cellManager = cellManager;
        this.bondingInfo = bondingInfo;
        this.isPureAtoms = isPureAtoms;
        this.box = box;
        this.space = box.getSpace();
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
            if (bondingInfo.skipBondedPair(isPureAtoms, atom1, atom2)) {
                continue;
            }
            Vector rij = space.makeVector();
            rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
            consumer.accept(atom2, rij);
        }

        int iCell = atomCell[iAtom];
        for (int cellOffset : cellOffsets) {
            int jCell = iCell + cellOffset;
            Vector jbo = boxOffsets[jCell];
            jCell = wrapMap[jCell];
            for (int j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                IAtom atom2 = atoms.get(j);
                if (bondingInfo.skipBondedPair(isPureAtoms, atom1, atom2)) {
                    continue;
                }

                Vector rij = space.makeVector();
                rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                rij.PE(jbo);
                consumer.accept(atom2, rij);
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
            if (bondingInfo.skipBondedPair(isPureAtoms, atom1, atom2)) {
                continue;
            }
            Vector rij = space.makeVector();
            rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
            consumer.accept(atom2, rij);
        }

        for (int cellOffset : cellOffsets) {
            int jCell = iCell - cellOffset;
            Vector jbo = boxOffsets[jCell];
            jCell = wrapMap[jCell];
            for (int j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                IAtom atom2 = atoms.get(j);
                if (bondingInfo.skipBondedPair(isPureAtoms, atom1, atom2)) {
                    continue;
                }

                Vector rij = space.makeVector();
                rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                rij.PE(jbo);
                consumer.accept(atom2, rij);
            }
        }
    }

    public interface SuperNbrConsumer {
        double accept(IAtom atom1, IAtom atom2, Vector rij);
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
            if (j == iAtom) { continue; }
            IAtom atom2 = atoms.get(j);
            if (bondingInfo.skipBondedPair(isPureAtoms, atom1, atom2)) {
                continue;
            }
            Vector rij = space.makeVector();
            rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
            sum += consumer.accept(atom1, atom2, rij);
        }

        for (int cellOffset : allCellOffsets) {
            int jCell = iCell + cellOffset;
            Vector jbo = boxOffsets[jCell];
            jCell = wrapMap[jCell];
            for (int j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                IAtom atom2 = atoms.get(j);
                if (bondingInfo.skipBondedPair(isPureAtoms, atom1, atom2)) {
                    continue;
                }

                Vector rij = space.makeVector();
                rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                rij.PE(jbo);
                sum += consumer.accept(atom1, atom2, rij);
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
//        int[] cellOffsets = cellManager.getCellOffsets();
        int[] allCellOffsets = cellManager.allCellOffsets;
        int[] wrapMap = cellManager.getWrapMap();
        int[] cellLastAtom = cellManager.getCellLastAtom();
        IAtom atom1 = atoms.get(iAtom);
        int iCell = atomCell[iAtom];

        for (int j = cellLastAtom[iCell]; j > -1; j = cellNextAtom[j]) {
            if (j == iAtom) { continue; }
            IAtom atom2 = atoms.get(j);
            if (bondingInfo.skipBondedPair(isPureAtoms, atom1, atom2)) {
                continue;
            }
            Vector rij = space.makeVector();
            rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
            consumer.accept(atom2, rij);
        }

        for (int cellOffset : allCellOffsets) {
            int jCell = iCell + cellOffset;
            Vector jbo = boxOffsets[jCell];
            jCell = wrapMap[jCell];
            for (int j = cellLastAtom[jCell]; j > -1; j = cellNextAtom[j]) {
                IAtom atom2 = atoms.get(j);
                if (bondingInfo.skipBondedPair(isPureAtoms, atom1, atom2)) {
                    continue;
                }

                Vector rij = space.makeVector();
                rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                rij.PE(jbo);
                consumer.accept(atom2, rij);
            }
        }
    }

}
