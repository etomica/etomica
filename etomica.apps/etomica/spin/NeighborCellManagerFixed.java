package etomica.spin;

import etomica.Atom;
import etomica.Phase;
import etomica.PhaseCellManager;
import etomica.PhaseEvent;
import etomica.PhaseListener;
import etomica.SimulationEvent;
import etomica.Space;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.lattice.CellLattice;
import etomica.lattice.RectangularLattice;
import etomica.nbr.site.AtomSite;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on May 23, 2005 by kofke
 */
public class NeighborCellManagerFixed implements PhaseCellManager {

    private final CellLattice lattice;
    private final Space space;
    private final AtomIteratorListSimple atomIterator;
    private final RectangularLattice.Iterator siteIterator;
    
    /**
     * Constructs manager for neighbor cells in the given phase.  The number of
     * cells in each dimension is given by nCells. Position definition for each
     * atom is that given by its type (it is set to null in this class).
     */
    public NeighborCellManagerFixed(Phase phase, int nCells) {
        space = phase.space();
        atomIterator = new AtomIteratorListSimple(phase.speciesMaster().atomList);

        lattice = new CellLattice(phase.boundary().dimensions(), AtomSite.FACTORY);
        int[] size = new int[space.D()];
        for(int i=0; i<space.D(); i++) size[i] = nCells;
        lattice.setSize(size);
        siteIterator = new RectangularLattice.Iterator(space.D());
        siteIterator.setLattice(lattice);

        //listener to phase to detect addition of new SpeciesAgent
        //or new atom
        phase.speciesMaster.addListener(new PhaseListener() {
            public void actionPerformed(SimulationEvent evt) {
                actionPerformed((PhaseEvent)evt);
            }
            public void actionPerformed(PhaseEvent evt) {
                if(evt.type() == PhaseEvent.ATOM_ADDED) {
                    Atom atom = evt.atom();
                    //new species agent requires another list in each cell
                    assignCell(atom);
                }
            }
        });
    }

    public CellLattice getLattice() {
        return lattice;
    }
    
    /**
     * Assigns cells to all leaf atoms in the phase.
     */
    public void assignCellAll() {
        atomIterator.reset();
        siteIterator.reset();
        while(atomIterator.hasNext()) {
            ((AtomSite)siteIterator.next()).setAtom(atomIterator.nextAtom());
        }
    }

    public void assignCell(Atom atom) {
        // TODO Auto-generated method stub

    }

}
