package etomica.nbr.cell.molecule;

import etomica.lattice.CellLattice;

/**
 * Interface for moleculeset iterators that use cell-based neighbor
 * lists for iteration.  Defines a method that returns the iterator
 * that loops over the cells neighboring a given cell.  Access to
 * the cell iterator is needed to adjust the neighbor range defining
 * the neighbor cells.
 *
 * @author Tai Boon Tan
 *
 */
public interface MoleculesetIteratorCellular {

    /**
     * @return Returns the neighborCell Iterator.
     */
    public CellLattice.NeighborIterator getNbrCellIterator();
}