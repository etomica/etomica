/*
 * Created on Jan 7, 2005
 */
package etomica.nbr.cell;

import etomica.lattice.CellLattice;

/**
 * Interface for atomset iterators that use cell-based neighbor
 * lists for iteration.  Defines a method that returns the iterator
 * that loops over the cells neighboring a given cell.  Access to
 * the cell iterator is needed to adjust the neighbor range defining
 * the neighbor cells.
 */
public interface AtomsetIteratorCellular {

    /**
     * @return Returns the neighborCell Iterator.
     */
    public CellLattice.NeighborIterator getNbrCellIterator();
}