package etomica;

import etomica.lattice.CellLattice;

/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Apr 15, 2005 by kofke
 */
public interface PhaseCellManager {

    public CellLattice getLattice();

    /**
     * Assigns cells to all molecules in the phase.
     */
    public void assignCellAll();

    /**
     * Assigns the cell for the given atom.
     * @param atom
     */
    public void assignCell(Atom atom);
}