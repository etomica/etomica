/*
 * History
 * Created on Nov 21, 2004 by kofke
 */
package etomica.nbratom.cell;

import etomica.Atom;
import etomica.Phase;
import etomica.PhaseEvent;
import etomica.PhaseListener;
import etomica.SimulationEvent;
import etomica.Space;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomPositionDefinitionSimple;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.lattice.CellLattice;

/**
 * Class that defines and manages construction and use of lattice of cells 
 * for cell-based neighbor listing.
 */

//TODO modify assignCellAll to loop through cells to get all atoms to be assigned
//no need for index when assigning cell
//different iterator needed

public class NeighborCellManager {

    private final CellLattice lattice;
    private final Space space;
    private final AtomIteratorTree atomIterator;
    private final AtomPositionDefinition positionDefinition;
    
    /**
     * Constructs manager for neighbor cells in the given phase.  The number of
     * cells in each dimension is given by nCells. 
     */
    public NeighborCellManager(Phase phase, int nCells) {
        this(phase,nCells,new AtomPositionDefinitionSimple());
    }
    
    public NeighborCellManager(Phase phase, int nCells, AtomPositionDefinition positionDefinition) {
        this.positionDefinition = positionDefinition;
        space = phase.space();
        atomIterator = new AtomIteratorTree();
        atomIterator.setDoAllNodes(true);
        atomIterator.setRoot(phase.speciesMaster);

        lattice = new CellLattice(phase.boundary().dimensions(), NeighborCell.FACTORY);
        phase.setLattice(lattice);
        int[] size = new int[space.D()];
        for(int i=0; i<space.D(); i++) size[i] = nCells;
        lattice.setSize(size);

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
                    if(atom.type.getNbrManagerAgent().getPotentials().length > 0) {
                        assignCell(atom);
                    }
                }
            }
        });
    }

    /**
     * Assigns cells to all molecules in the phase.
     */
    public void assignCellAll() {
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom atom = atomIterator.nextAtom();
            if (atom.type.getNbrManagerAgent().getPotentials().length > 0) {
                assignCell(atom);
            }
        }
    }
    
    /**
     * Assigns the cell for the given atom.
     * @param atom
     */
    public void assignCell(Atom atom) {
        AtomSequencerCell seq = (AtomSequencerCell)atom.seq;
        NeighborCell newCell = (NeighborCell)lattice.site(positionDefinition.position(atom));
        if(newCell != seq.cell) {assignCell(seq, newCell);}
    }
    
    /**
     * Assigns atom sequencer to given cell in the list of the given index.
     */
    public void assignCell(AtomSequencerCell seq, NeighborCell newCell) {
        if(seq.cell != null) seq.cell.occupants().remove(seq.nbrLink);
        seq.cell = newCell;
        if(newCell != null) {
            newCell.occupants().add(seq.nbrLink);
        }
    }//end of assignCell
    
}//end of NeighborCellManager
