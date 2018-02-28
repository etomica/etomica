/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.cell;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomSetSinglet;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.box.Box;
import etomica.box.BoxCellManager;
import etomica.integrator.mcmove.MCMove;
import etomica.integrator.mcmove.MCMoveEvent;
import etomica.integrator.mcmove.MCMoveTrialCompletedEvent;
import etomica.integrator.mcmove.MCMoveTrialFailedEvent;
import etomica.lattice.CellLattice;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.simulation.Simulation;
import etomica.space.*;
import etomica.util.Debug;
import etomica.util.IListener;

/**
 * Class that defines and manages construction and use of lattice of cells 
 * for cell-based neighbor listing.
 *
 * Before starting the simulation, it is necessary call makeMCMoveListener and
 * register the result with the MCMoveEventManager of the Integrator.
 */

//TODO modify assignCellAll to loop through cells to get all atoms to be assigned
//no need for index when assigning cell
//different iterator needed

public class NeighborCellManager implements BoxCellManager, BoundaryEventListener, AtomLeafAgentManager.AgentSource<Cell> {

    protected final Simulation sim;
    protected final CellLattice lattice;
    protected final IMoleculePositionDefinition positionDefinition;
    protected final Box box;
    protected int cellRange = 2;
    protected double range;
    protected final AtomLeafAgentManager<Cell> agentManager;
    protected boolean doApplyPBC;
    protected final Vector v;
    protected final int[] numCells;
    protected boolean suppressBoxLengthWarning;
    
    /**
     * Constructs manager for neighbor cells in the given box.  The number of
     * cells in each dimension is given by nCells. Position definition for each
     * atom is that given by its type (it is set to null in this class).
     */
    public NeighborCellManager(Simulation sim, Box box, double potentialRange, Space _space) {
        this(sim, box, potentialRange, null, _space);
    }
    
    /**
     * Construct manager for neighbor cells in the given box.  The number
     * of cells in each dimension is given by nCells.  Position definition is
     * used to determine the cell a given atom is in; if null, the position
     * definition given by the atom's type is used.  Position definition is
     * declared final.
     */
    public NeighborCellManager(Simulation sim, Box box, double potentialRange, IMoleculePositionDefinition positionDefinition, Space space) {
        this.positionDefinition = positionDefinition;
        this.box = box;
        this.sim = sim;
        numCells = new int[space.D()];

        lattice = new CellLattice(space, box.getBoundary().getBoxSize(), Cell.FACTORY);
        setPotentialRange(potentialRange);
        v = space.makeVector();
        agentManager = new AtomLeafAgentManager<Cell>(this,box);
        doApplyPBC = false;
    }

    /**
     * @param doSuppress warnings about box length being too small.  This is
     *                   useful if the potential range is short enough, but the
     *                   neighbor listing padding does not fit in the box.
     */
    public void setSuppressBoxLengthWarning(boolean doSuppress) {
        suppressBoxLengthWarning = doSuppress;
    }
    
    public void setDoApplyPBC(boolean newDoApplyPBC) {
        doApplyPBC = newDoApplyPBC;
    }
    
    public boolean getDoApplyPBC() {
        return doApplyPBC;
    }

    public CellLattice getLattice() {
        return lattice;
    }

    /**
     * Sets the potential range to the given value.  Cells are made large 
     * enough so that {@code cellRange*cellSize > potentialRange}.
     */
    public void setPotentialRange(double newRange) {
        range = newRange;
        if (checkDimensions() && agentManager != null) {
            assignCellAll();
            lattice.setNeighborRange(newRange);
        }
    }
    
    /**
     * Returns the potential range.
     */
    public double getPotentialRange() {
        return range;
    }
    
    /**
     * Returns the cellRange.
     */
    public int getCellRange() {
        return cellRange;
    }

    /**
     * Sets the cell range to the given value.  Cells are made large 
     * enough so that cellRange*cellSize > potentialRange
     */
    public void setCellRange(int newCellRange) {
        cellRange = newCellRange;
        checkDimensions();
    }

    public void boundaryInflate(BoundaryEvent e) {
        checkDimensions();
        // we need to reassign cells even if checkDimensions didn't resize
        // the lattice.  If the box size changed, the cell size changed,
        // and the atom assignments need to change too.
        //FIXME but only if we have multi-atomic molecules.  For monatomic
        // molecules, we would only need to call this if the lattice size
        // changes
        boolean savedDoApplyPBC = doApplyPBC;
        doApplyPBC = true;
        assignCellAll();
        doApplyPBC = savedDoApplyPBC;
    }
    
    /**
     * Checks the box's dimensions to make sure the number of cells is 
     * appropriate.
     */
    protected boolean checkDimensions() {
        if (range == 0) {
            // simulation is still being constructed, don't try to do anything useful
            return false;
        }
        Vector dimensions = box.getBoundary().getBoxSize();
        lattice.setDimensions(dimensions);
        int[] oldSize = lattice.getSize();
        boolean latticeNeedsUpdate = false;
        for (int i=0; i<numCells.length; i++) {
            numCells[i] = (int)Math.floor(cellRange*dimensions.getX(i)/range);
            if (numCells[i] > 100) {
                // too many cells will cause us to run out of memory.
                // >100 cells is generally not useful
                // we usually end up here because we will later increase the
                // potential range
                if (Debug.ON) System.err.println("capping "+i+" cells at 100");
                numCells[i] = 100;
            }
            else if (numCells[i] < cellRange*2+1) {
                // the box is too small for the lattice, but things might still be OK.
                // it's possible the box is big enough for the potential range, but
                // the extra padding that's needed for happy cells is missing.
                // capping the number of cells means we can proceed and just don't
                // get any advantage from using cells in this direction.  ideally
                // we would make the neighbor iterator non-wrapping in this direction
                // and use 1 cell.
                if (Debug.ON) System.err.println("bumping number of cells in direction "+i+" from "+numCells[i]+" to "+(cellRange*2+1));
                numCells[i] = cellRange*2+1;
                if (range > dimensions.getX(i)/2 && !suppressBoxLengthWarning) {
                    // box was too small for the potentials too.  doh.
                    // Perhaps the direction is not periodic or we're in the middle
                    // of multiple changes which will (in the end) be happy.
                    System.err.println("range is greater than half the box length in direction "+i);
                }
            }
            latticeNeedsUpdate = latticeNeedsUpdate || oldSize[i] != numCells[i];
        }

        //only update the lattice (expensive) if the number of cells changed
        if (latticeNeedsUpdate) {
            lattice.setSize(numCells);
            return true;
        }
        return false;
    }
    
    /**
     * Assigns cells to all interacting atoms in the box.  Interacting atoms
     * are those that have one or more potentials that act on them.  
     */
    public void assignCellAll() {
        // ensure that any changes to cellRange, potentialRange and boundary
        // dimension take effect.  checkDimensions can call us, but if that
        // happens, our call into checkDimensions should 
        checkDimensions();

        Object[] allCells = lattice.sites();
        for (int i=0; i<allCells.length; i++) {
            ((Cell)allCells[i]).occupants().clear();
        }

        IAtomList leafList = box.getLeafList();
        int count = leafList.size();
        for (int i=0; i<count; i++) {
            IAtom atom = leafList.get(i);
            assignCell(atom);
        }
    }
    
    public Cell getCell(IAtom atom) {
        return agentManager.getAgent(atom);
    }

    /**
     * Assigns the cell for the given atom.  The atom will be listed in the
     * cell's atom list and the cell with be associated with the atom via
     * agentManager.
     */
    public void assignCell(IAtom atom) {
        Cell atomCell;
        if (doApplyPBC) {
            v.E(atom.getPosition());
            v.PE(box.getBoundary().centralImage(v));
            atomCell = (Cell)lattice.site(v);
        }
        else {
            atomCell = (Cell)lattice.site(atom.getPosition());
        }
        atomCell.addAtom(atom);
        agentManager.setAgent(atom, atomCell);
    }

    public IListener<MCMoveEvent> makeMCMoveListener() {
        return new MyMCMoveListener(box,this);
    }

    /**
     * Returns the cell containing the given atom.  The atom is added to the
     * cell's atom list.
     */
    public Cell makeAgent(IAtom atom, Box agentBox) {
        // if we have no cells, there's no point in trying here.  cell assignment will happen later
        if (numCells[0] == 0) return null;
        Vector position = atom.getPosition();
        v.E(position);
        v.PE(box.getBoundary().centralImage(position));
        Cell atomCell = (Cell)lattice.site(v);
        atomCell.addAtom(atom);
        if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(atom))) {
            System.out.println("assigning new "+atom+" at "+position+" to "+atomCell);
        }
        return atomCell;
    }

    /**
     * Removes the given atom from the cell.
     */
    public void releaseAgent(Cell cell, IAtom atom, Box agentBox) {
        cell.removeAtom(atom);
    }

    private static class MyMCMoveListener implements IListener<MCMoveEvent>, java.io.Serializable {
        public MyMCMoveListener(Box box, NeighborCellManager manager) {
            this.box = box;
            neighborCellManager = manager;
        }

        public void actionPerformed(MCMoveEvent evt) {
            if (evt instanceof MCMoveTrialCompletedEvent && ((MCMoveTrialCompletedEvent)evt).isAccepted()) {
                return;
            }
            if (evt instanceof MCMoveTrialFailedEvent) {
                return;
            }

            MCMove move = evt.getMCMove();
            AtomIterator iterator = move.affectedAtoms(box);
            iterator.reset();
            for (IAtom atom = iterator.nextAtom(); atom != null; atom = iterator.nextAtom()) {
                updateCell(atom);
            }
        }

        private void updateCell(IAtom atom) {
            Boundary boundary = box.getBoundary();
            Cell cell = neighborCellManager.getCell(atom);
            cell.removeAtom(atom);
            atom.getPosition().PE(boundary.centralImage(atom.getPosition()));
            neighborCellManager.assignCell(atom);
        }
        
        private static final long serialVersionUID = 1L;
        private final Box box;
        private final NeighborCellManager neighborCellManager;
    }
}
