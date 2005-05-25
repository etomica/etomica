/*
 * History
 * Created on Nov 21, 2004 by kofke
 */
package etomica.nbratom.cell;

import etomica.Atom;
import etomica.AtomIterator;
import etomica.AtomTypeLeaf;
import etomica.Phase;
import etomica.PhaseCellManager;
import etomica.PhaseEvent;
import etomica.PhaseListener;
import etomica.SimulationEvent;
import etomica.Space;
import etomica.action.AtomActionTranslateBy;
import etomica.action.AtomGroupAction;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.data.DataSourceCOM;
import etomica.integrator.MCMove;
import etomica.integrator.mcmove.MCMoveEvent;
import etomica.integrator.mcmove.MCMoveListener;
import etomica.lattice.CellLattice;
import etomica.space.Boundary;
import etomica.space.Vector;

/**
 * Class that defines and manages construction and use of lattice of cells 
 * for cell-based neighbor listing.
 */

//TODO modify assignCellAll to loop through cells to get all atoms to be assigned
//no need for index when assigning cell
//different iterator needed

public class NeighborCellManager implements PhaseCellManager {

    private final CellLattice lattice;
    private final Space space;
    private final AtomIteratorTree atomIterator;
    private final AtomPositionDefinition positionDefinition;
    private final Phase phase;

    
    /**
     * Constructs manager for neighbor cells in the given phase.  The number of
     * cells in each dimension is given by nCells. Position definition for each
     * atom is that given by its type (it is set to null in this class).
     */
    public NeighborCellManager(Phase phase, int nCells) {
        this(phase,nCells, null);
    }
    
    /**
     * Construct manager for neighbor cells in the given phase.  The number
     * of cells in each dimension is given by nCells.  Position definition is
     * used to determine the cell a given atom is in; if null, the position
     * definition given by the atom's type is used.  Position definition is
     * declared final.
     */
    public NeighborCellManager(Phase phase, int nCells, AtomPositionDefinition positionDefinition) {
        this.positionDefinition = positionDefinition;
        this.phase = phase;
        space = phase.space();
        atomIterator = new AtomIteratorTree();
        atomIterator.setDoAllNodes(true);
        atomIterator.setRoot(phase.speciesMaster);

        lattice = new CellLattice(phase.boundary().dimensions(), NeighborCell.FACTORY);
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

    public CellLattice getLattice() {
        return lattice;
    }
    
    /**
     * Assigns cells to all molecules in the phase.
     */
    public void assignCellAll() {
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom atom = atomIterator.nextAtom();
            if (atom.type.isInteracting()  && (atom.type instanceof AtomTypeLeaf && ((AtomTypeLeaf)atom.type).getMass()!=Double.POSITIVE_INFINITY ||
                    ((AtomSequencerCell)atom.seq).cell == null)) {
                assignCell(atom);
            }
        }
    }
    
    /**
     * Assigns the cell for the given atom.
     * @param atom
     */
    public void assignCell(Atom atom) {
        Vector position = (positionDefinition != null) ?
                positionDefinition.position(atom) :
                    atom.type.getPositionDefinition().position(atom);
        ((NeighborCell)lattice.site(position)).addAtom(atom);
    }
    
    public MCMoveListener makeMCMoveListener() {
        return new MyMCMoveListener();
    }

    
    private class MyMCMoveListener implements MCMoveListener {
        public MyMCMoveListener() {
            treeIterator = new AtomIteratorTree();
            treeIterator.setDoAllNodes(true);
            moleculePosition = new DataSourceCOM(space);
            translator = new AtomActionTranslateBy(space);
            moleculeTranslator = new AtomGroupAction(translator);
        }
        
        public void actionPerformed(SimulationEvent evt) {
            actionPerformed((MCMoveEvent)evt);
        }

        public void actionPerformed(MCMoveEvent evt) {
            if (!evt.isTrialNotify && evt.wasAccepted) {
                return;
            }
            MCMove move = evt.mcMove;
            AtomIterator iterator = move.affectedAtoms(phase);
            iterator.reset();
            while (iterator.hasNext()) {
                Atom atom = iterator.nextAtom();
                if (!atom.node.isLeaf()) {
                    treeIterator.setRoot(atom);
                    treeIterator.reset();
                    while (treeIterator.hasNext()) {
                        Atom childAtom = treeIterator.nextAtom();
                        updateCell(childAtom);
                    }
                }
                else {
                    updateCell(atom);
                }
            }
        }

        private void updateCell(Atom atom) {
            Boundary boundary = phase.boundary();
            if (((AtomSequencerCell)atom.seq).cell != null) {
                if (!atom.node.isLeaf()) {
                    Vector shift = boundary.centralImage(moleculePosition.position(atom));
                    if (!shift.isZero()) {
                        translator.setTranslationVector(shift);
                        moleculeTranslator.actionPerformed(atom);
                    }
                }
                else {
                    Vector shift = boundary.centralImage(atom.coord.position());
                    if (!shift.isZero()) {
                        atom.coord.position().PE(shift);
                    }
                }
                assignCell(atom);
            }
        }
        
        private final AtomIteratorTree treeIterator;
        private final AtomPositionDefinition moleculePosition;
        private final AtomActionTranslateBy translator;
        private final AtomGroupAction moleculeTranslator;
    }
}//end of NeighborCellManager
