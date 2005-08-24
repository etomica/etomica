package etomica.spin;

import etomica.Phase;
import etomica.PhaseCellManager;
import etomica.PhaseEvent;
import etomica.PhaseListener;
import etomica.SimulationEvent;
import etomica.atom.Atom;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.AtomIteratorListTabbed;
import etomica.lattice.CellLattice;
import etomica.lattice.RectangularLattice;
import etomica.nbr.site.AtomSite;
import etomica.space.Space;

/**
 * Neighbor manager for system in which there is an unchanging, one-to-one
 * correspondence between atoms and lattice cells. Each leaf atom is assigned to
 * its own cell at the beginning of the simulation, and the cell-atom assignment
 * doesn't change through the course of the simulation. Atom neighbors are
 * determined by the structure and neighbor definition associated with the cell
 * lattice, and do not change through the simulation.
 * <p>
 * The expected use for this neighbor structure is in modeling of
 * lattice systems such as the Ising model in which the coordinate represents
 * something other than the spatial position of the atom; can also be applied to
 * diffusionless crystalline solids, such as the valence-force field model,
 * single-occupancy crystals, and perhaps high-density unconstrained crystals.
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on May 23, 2005 by kofke
 */
public class NeighborCellManagerFixed implements PhaseCellManager {

    /**
     * Constructs manager for neighbor cells in the given phase. The number of
     * cells in each dimension is given by nCells. Position definition for each
     * atom is that given by its type (it is set to null in this class).
     */
    public NeighborCellManagerFixed(final Phase phase, int nCells) {
        space = phase.space();
        atomIterator = new AtomIteratorListTabbed(
                phase.getSpeciesMaster().atomList);

        lattice = new CellLattice(phase.boundary().dimensions(),
                AtomSite.FACTORY);
        int[] size = new int[space.D()];
        for (int i = 0; i < space.D(); i++)
            size[i] = nCells;
        lattice.setSize(size);
        siteIterator = new RectangularLattice.Iterator(space.D());
        siteIterator.setLattice(lattice);

        //listener to phase to detect addition of new SpeciesAgent
        //or new atom
        ((SpeciesRoot)phase.getSpeciesMaster().node.parentGroup()).addListener(new PhaseListener() {

            public void actionPerformed(SimulationEvent evt) {
                if (((PhaseEvent)evt).phase() == phase) {
                    actionPerformed((PhaseEvent)evt);
                }
            }

            public void actionPerformed(PhaseEvent evt) {
                if (evt.type() == PhaseEvent.ATOM_ADDED) {
//                    throw new RuntimeException(
//                            "Neighbor cell manager cannot handle addition of atoms");
                }
            }
        });
    }

    /**
     * Returns the lattice used to define the neighbor structure.
     */
    public CellLattice getLattice() {
        return lattice;
    }

    /**
     * Assigns cells to all leaf atoms in the phase.
     */
    public void assignCellAll() {
        atomIterator.reset();
        siteIterator.reset();
        while (atomIterator.hasNext()) {
            ((AtomSite) siteIterator.next()).setAtom(atomIterator.nextAtom());
        }
    }

    /**
     * Should not be called, because cell assignments can be done only once, and
     * is handled by assignCellAll.
     * 
     * @throws RuntimeException
     *             if invoked
     */
    public void assignCell(Atom atom) {
        throw new RuntimeException(
                "Cell assignments can be made only via assignCellAll, at the beginning of the simulation");
    }

    private final CellLattice lattice;
    private final Space space;
    private final AtomIteratorListTabbed atomIterator;
    private final RectangularLattice.Iterator siteIterator;

}
