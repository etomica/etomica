package etomica.nbr.cell;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.lattice.CellLattice;

import java.util.function.BiConsumer;

public class ApiAACellFast {
    private final Box box;
    private final CellLattice lattice;

    public ApiAACellFast(Box box, CellLattice lattice) {
        this.box = box;
        this.lattice = lattice;
    }

    public void forEachPair(AtomPairConsumer action) {
        Cell[] sites = (Cell[]) lattice.sites();
        int[][] nbrCells = lattice.getUpNeighbors();
        for (int centralCellIdx = 0; centralCellIdx < sites.length; centralCellIdx++) {

            Cell centralCell = sites[centralCellIdx];
            IAtomList centralCellAtoms = centralCell.occupants();

            // loop over pairs within the cell
            for (int i = 0; i < centralCellAtoms.size(); i++) {
                for (int j = i + 1; j < centralCellAtoms.size(); j++) {
                    action.accept(centralCellAtoms.get(i), centralCellAtoms.get(j));
                }
            }

            for (int nbrCellIdx : nbrCells[centralCellIdx]) {
                Cell nbrCell = sites[nbrCellIdx];
                IAtomList nbrCellAtoms = nbrCell.occupants();

                for (int i = 0; i < centralCellAtoms.size(); i++) {
                    for (int j = 0; j < nbrCellAtoms.size(); j++) {
                        action.accept(centralCellAtoms.get(i), nbrCellAtoms.get(j));
                    }
                }

            }
        }
    }

    @FunctionalInterface
    public interface AtomPairConsumer extends BiConsumer<IAtom, IAtom> {
        @Override
        void accept(IAtom atom1, IAtom atom2);
    }
}
