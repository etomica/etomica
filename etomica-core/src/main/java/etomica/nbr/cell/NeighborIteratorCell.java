package etomica.nbr.cell;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.lattice.CellLattice;
import etomica.nbr.NeighborIterator;
import etomica.potential.IteratorDirective;

public class NeighborIteratorCell implements NeighborIterator {

    private final NeighborCellManager neighborCellManager;
    private final CellLattice lattice;

    public NeighborIteratorCell(NeighborCellManager neighborCellManager) {
        this.neighborCellManager = neighborCellManager;
        this.lattice = neighborCellManager.getLattice();
    }

    @Override
    public void forEachNeighbor(IAtom targetAtom, IteratorDirective.Direction direction, AtomPairConsumer upAction, AtomPairConsumer downAction) {

        Cell targetCell = neighborCellManager.getCell(targetAtom);
        Object[] sites = lattice.sites();
        IAtomList cellAtoms = targetCell.occupants();
        int[] upNbrCells = lattice.getUpNeighbors()[targetCell.getLatticeArrayIndex()];
        int[] downNbrCells = lattice.getDownNeighbors()[targetCell.getLatticeArrayIndex()];
        boolean doUp = direction != IteratorDirective.Direction.DOWN;
        boolean doDown = direction != IteratorDirective.Direction.UP;
        boolean upListNow = false;

        for (int i = 0; i < cellAtoms.size(); i++) {
            IAtom otherAtom = cellAtoms.get(i);
            if(otherAtom == targetAtom) {
                upListNow = true;
            } else {
                if(upListNow && doUp) {
                    upAction.accept(targetAtom, otherAtom);
                } else if (doDown) {
                    downAction.accept(otherAtom, targetAtom);
                }
            }
        }

        if (doUp) {
            for (int nbrIdx : upNbrCells) {
                Cell nbrCell = (Cell) sites[nbrIdx];
                for (int i = 0; i < nbrCell.occupants().size(); i++) {
                    upAction.accept(targetAtom, nbrCell.occupants().get(i));
                }
            }
        }

        if (doDown) {
            for (int nbrIdx : downNbrCells) {
                Cell nbrCell = (Cell) sites[nbrIdx];
                for (int i = 0; i < nbrCell.occupants().size(); i++) {
                    downAction.accept(nbrCell.occupants().get(i), targetAtom);
                }
            }
        }
    }
}
