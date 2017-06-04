/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.cell.molecule;

import etomica.atom.*;
import etomica.space.Boundary;
import etomica.space.BoundaryEvent;
import etomica.api.IBoundaryListener;
import etomica.box.Box;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.atom.MoleculePositionGeometricCenter;
import etomica.atom.IMoleculePositionDefinition;
import etomica.atom.iterator.AtomIterator;
import etomica.box.BoxCellManager;
import etomica.integrator.mcmove.MCMove;
import etomica.integrator.mcmove.MCMoveEvent;
import etomica.integrator.mcmove.MCMoveTrialCompletedEvent;
import etomica.lattice.CellLattice;
import etomica.space.Space;
import etomica.util.Debug;
import etomica.util.IEvent;
import etomica.util.IListener;

/**
 * Class that defines and manages construction and use of lattice of cells 
 * for cell-based neighbor listing.
 * 
 * 
 * @author Tai Boon Tan
 *
 */

//TODO modify assignCellAll to loop through cells to get all atoms to be assigned
//no need for index when assigning cell
//different iterator needed

public class NeighborCellManagerMolecular implements BoxCellManager, IBoundaryListener, MoleculeAgentManager.MoleculeAgentSource, java.io.Serializable {

    private static final long serialVersionUID = 1L;
    protected final Simulation sim;
    protected final CellLattice lattice;
    protected final IMoleculePositionDefinition positionDefinition;
    protected final Box box;
    protected int cellRange = 2;
    protected double range;
    protected final MoleculeAgentManager agentManager;
    protected boolean doApplyPBC;
    protected final Vector v;
    protected final int[] numCells;
    protected IMoleculePositionDefinition moleculeSite;
    protected Space space;
    
    /**
     * Constructs manager for neighbor cells in the given box.  The number of
     * cells in each dimension is given by nCells. Position definition for each
     * atom is that given by its type (it is set to null in this class).
     */
    public NeighborCellManagerMolecular(Simulation sim, Box box, double potentialRange, Space _space) {
        this(sim, box, potentialRange, null, _space);
    }
    
    /**
     * Construct manager for neighbor cells in the given box.  The number
     * of cells in each dimension is given by nCells.  Position definition is
     * used to determine the cell a given atom is in; if null, the position
     * definition given by the atom's type is used.  Position definition is
     * declared final.
     */
    public NeighborCellManagerMolecular(Simulation sim, Box box, double potentialRange, IMoleculePositionDefinition positionDefinition, Space space) {
        this.positionDefinition = positionDefinition;
        this.box = box;
        this.sim = sim;
        this.space = space;
        numCells = new int[space.D()];

        lattice = new CellLattice(space, box.getBoundary().getBoxSize(), CellMolecular.FACTORY);
        setPotentialRange(potentialRange);
        v = space.makeVector();
        agentManager = new MoleculeAgentManager(sim, box, this);
        doApplyPBC = false;
        moleculeSite = new MoleculePositionGeometricCenter(space);
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
     * enough so that cellRange*cellSize > potentialRange.
     */
    public void setPotentialRange(double newRange) {
        range = newRange;
        if (checkDimensions() && agentManager != null) {
            assignCellAll();
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
        assignCellAll();
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
                if (range > dimensions.getX(i)/2) {
                    // box was too small for the potentials too.  doh.
                    // Perhaps the direction is not periodic or we're in the middle
                    // of multiple changes which will (in the end) be happy.
                    System.err.println("<NeighborCellManagerMolecular> range is greater than half the box length in direction "+i + "; " + dimensions.getX(i)/2);
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
            ((CellMolecular)allCells[i]).occupants().clear();
        }
        
        IMoleculeList moleculeList = box.getMoleculeList();
        int count = moleculeList.getMoleculeCount();
        for (int i=0; i<count; i++) {
            IMolecule molecule = moleculeList.getMolecule(i);
            assignCell(molecule);
        }
    }
    
    public CellMolecular getCell(IMolecule molecule) {
        return (CellMolecular)agentManager.getAgent(molecule);
    }

    /**
     * Assigns the cell for the given atom.  The atom will be listed in the
     * cell's atom list and the cell with be associated with the atom via
     * agentManager.
     */
    public void assignCell(IMolecule molecule) {
        CellMolecular moleculeCell;
        if (doApplyPBC) {
            v.E(moleculeSite.position(molecule));
            v.PE(box.getBoundary().centralImage(v));
            moleculeCell = (CellMolecular)lattice.site(v);
        }
        else {
            moleculeCell = (CellMolecular)lattice.site(moleculeSite.position(molecule));
        }
        moleculeCell.addMolecule(molecule);
        agentManager.setAgent(molecule, moleculeCell);
    }
    
    public IListener makeMCMoveListener() {
        return new MyMCMoveListener(box,this, space);
    }
    
    public Class getMoleculeAgentClass() {
        return CellMolecular.class;
    }

    /**
     * Returns the cell containing the given atom.  The atom is added to the
     * cell's atom list.
     */
    public Object makeAgent(IMolecule molecule) {
        Vector position = moleculeSite.position(molecule);
        v.E(position);
        if (doApplyPBC) {
            v.PE(box.getBoundary().centralImage(position));
        }
        CellMolecular moleculeCell = (CellMolecular)lattice.site(v);
        moleculeCell.addMolecule(molecule);
//        if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(atom))) {
//            System.out.println("assigning new "+atom+" at "+position+" to "+atomCell);
//        }
        return moleculeCell;
    }

    /**
     * Removes the given atom from the cell.
     */
    public void releaseAgent(Object cell, IMolecule molecule) {
        ((CellMolecular)cell).removeMolecule(molecule);
    }
    
    private static class MyMCMoveListener implements IListener, java.io.Serializable {
        public MyMCMoveListener(Box box, NeighborCellManagerMolecular manager, Space space) {
            this.box = box;
            neighborCellManager = manager;
            moleculePosition = new MoleculePositionGeometricCenter(space);
        }
/*        
 * 
 */
        public void actionPerformed(IEvent evt) {
        	
        	
            if (evt instanceof MCMoveTrialCompletedEvent && ((MCMoveTrialCompletedEvent)evt).isAccepted()) {
            	
                return;
            }

            if (evt instanceof MCMoveEvent) {
                MCMove move = ((MCMoveEvent)evt).getMCMove();
                AtomIterator iterator = move.affectedAtoms(box);
                iterator.reset();
                if(iterator == null){
                	throw new RuntimeException("<NeighborCellManagerMolecular> MoleculeIterator!!!!!!!!!!!!");
                }
                for (IAtom atom = iterator.nextAtom(); atom != null; atom = iterator.nextAtom()) {
                	if(atom.getIndex() == 1);
                    updateCell(atom.getParentGroup());
                }
            }
        }

        private void updateCell(IMolecule molecule) {
            Boundary boundary = box.getBoundary();
            CellMolecular cell = neighborCellManager.getCell(molecule);
            cell.removeMolecule(molecule);
            boundary.nearestImage(moleculePosition.position(molecule));
            neighborCellManager.assignCell(molecule);
        }
        
        private static final long serialVersionUID = 1L;
        private final Box box;
        private final NeighborCellManagerMolecular neighborCellManager;
        private IMoleculePositionDefinition moleculePosition;
       
    }

	
}
