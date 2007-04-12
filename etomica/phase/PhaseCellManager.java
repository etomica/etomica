package etomica.phase;

import etomica.atom.IAtom;
import etomica.lattice.CellLattice;

/**
 * Interface for class that handles assignment of atoms to cells in a phase.
 * This facility is needed by neighbor-listing schemes, though it may find use for
 * other purposes.
 *
 * @author David Kofke and Andrew Schultz
 *
 */
public interface PhaseCellManager {

    /**
     * Returns the lattice that defines the cell arrangement.
     */
    public CellLattice getLattice();

    /**
     * Assigns cells to all molecules in the phase.
     */
    public void assignCellAll();

    /**
     * Assigns the cell for the given atom.
     */
    public void assignCell(IAtom atom);
}