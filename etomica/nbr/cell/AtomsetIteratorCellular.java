/*
 * Created on Jan 7, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.nbr.cell;

import etomica.lattice.CellLattice;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public interface AtomsetIteratorCellular {

    /**
     * @return Returns the cellIterator.
     */
    public CellLattice.NeighborIterator getNbrCellIterator();
}