/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.action.AtomsetDetect;
import etomica.atom.Atom;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.AtomToArrayListFixed;
import etomica.atom.iterator.ApiInnerFixed;
import etomica.atom.iterator.ApiSequence1A;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayList;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.atom.iterator.AtomPairIterator;
import etomica.atom.iterator.AtomsetIteratorPDT;
import etomica.atom.iterator.IteratorDirective;
import etomica.atom.iterator.IteratorDirective.Direction;
import etomica.lattice.CellLattice;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;
import etomica.space.BoundaryPeriodic;

/**
 * Generates pairs that are cell-based neighbors of a specific Atom. Iteration is
 * performed using cell lists, which defines the neighboring molecules.
 * Direction is related to ordering of cells and, within a cell, ordering of
 * molecules in cell's occupant list.
 */

public class Api1ACell implements AtomsetIteratorPDT, AtomsetIteratorCellular, 
        AtomPairIterator, java.io.Serializable {
    
    /**
     * Constructor makes iterator that must have phase specified and then be
     * reset() before iteration.
     * 
     * @param D
     *            the dimension of the space of the simulation (used to
     *            construct cell iterators)
     * @param range
     *            the distance within which pairs of atoms are considered
     *            neighbors. Used to define neighbor cells; some iterates may
     *            exceed this separation
     *  
     */
    public Api1ACell(int D, double range, PhaseAgentManager agentManager) {
        neighborIterator = new CellLattice.NeighborIterator(D, range);
        aiOuter = new AtomIteratorSinglet();
        aiSeq = new AtomIteratorArrayListSimple();
        //this iterator is used to loop through list of occupants of atoms's cell;
        //construct with AtomToLinker that gives appropriate linker
//        MyAtomToLinker atomToLinker = new MyAtomToLinker();
        atomToArrayListFixed = new AtomToArrayListFixed();
        aiSeqDirectableUp = new AtomIteratorArrayList(IteratorDirective.Direction.UP, 1, atomToArrayListFixed, atomToArrayListFixed);
        aiSeqDirectableDn = new AtomIteratorArrayList(IteratorDirective.Direction.DOWN, 1, atomToArrayListFixed, atomToArrayListFixed);
        nbrCellListIteratorInner = new ApiSequence1A(aiSeqDirectableUp,aiSeqDirectableDn); //used only by allAtoms
        nbrCellListIteratorUp = new ApiInnerFixed(aiOuter, aiSeq);//used only by allAtoms
        nbrCellListIteratorDn = new ApiInnerFixed(aiOuter, aiSeq, true);//used only by allAtoms
        latticeIndex = new int[D];

        neighborIterator.setDirection(null);
        phaseAgentManager = agentManager;
        setPhase(null);
	}

	public void setPhase(Phase phase) {
        if(phase != null) {
            NeighborCellManager[] cellManagers = (NeighborCellManager[])phaseAgentManager.getAgents();
            lattice = cellManagers[phase.getIndex()].getLattice();
            neighborIterator.setLattice(lattice);
            neighborIterator.setPeriod(phase.getBoundary().getDimensions());
            neighborIterator.setPeriodicity(((BoundaryPeriodic)phase.getBoundary()).getPeriodicity());
        }
        else {
            lattice = null;
        }
	}

    /**
     * Performs action on all iterates.
     */
    public void allAtoms(AtomsetAction action) {
        if(pair.atom0 == null) return;
        aiOuter.setAtom(targetAtom);
        neighborIterator.checkDimensions();
        lattice.latticeIndex(centralCell.latticeArrayIndex,latticeIndex);
        
        //get pairs in targetMolecule's cell
        atomToArrayListFixed.setArrayList(centralCell.occupants());
        nbrCellListIteratorInner.setAtom(targetAtom);
        nbrCellListIteratorInner.allAtoms(action);

        //loop over neighbor cells up from central cell
        neighborIterator.setSite(latticeIndex);
        if (direction != IteratorDirective.Direction.DOWN) {
            neighborIterator.setDirection(IteratorDirective.Direction.UP);
            neighborIterator.reset();
            while(neighborIterator.hasNext()) {
                Cell neighborCell = (Cell)neighborIterator.next();
                aiSeq.setList(neighborCell.occupants());
                if(neighborCell.occupants().size() > 0) nbrCellListIteratorUp.allAtoms(action);
            }
        }

        //loop over neighbor cells down from central cell
        if (direction != IteratorDirective.Direction.UP) {
            neighborIterator.setDirection(IteratorDirective.Direction.DOWN);
            neighborIterator.reset();
            while(neighborIterator.hasNext()) {
                Cell neighborCell = (Cell)neighborIterator.next(); 
                aiSeq.setList(neighborCell.occupants());
                if(neighborCell.occupants().size() > 0) nbrCellListIteratorDn.allAtoms(action);
            }
        }
    }//end of allAtoms
    
	/**
	 * Returns the number of atom pairs the iterator will return if
	 * reset and iterated in its present state.
	 */
	public int size() {
        AtomsetCount counter = new AtomsetCount();
        allAtoms(counter);
        return counter.callCount();
	}
	
	/**
	 * Indicates whether the given atom pair will be among the iterates
	 * given by the iterator if reset in its present state.  True only
	 * if an iterated pair would match the atoms as ordered in the given
	 * array.
	 */
	public boolean contains(AtomSet atoms) {
        if(!(atoms instanceof AtomPair) || ((AtomPair)atoms).atom0 == ((AtomPair)atoms).atom1) return false;
        AtomsetDetect detector = new AtomsetDetect(atoms);
        allAtoms(detector);
        return detector.detectedAtom();
	}

    public boolean hasNext() {
        return aiInner.hasNext();
    }
    
    public AtomSet next() {
        return nextPair();
    }
    
    public AtomPair nextPair() {
        if (!hasNext()) return null;
        if (upListNow) {
            pair.atom1 = aiInner.nextAtom();
            if (needUpdate) {
        		needUpdate = false;
    			pair.atom0 = targetAtom;
            }
        }
        else {
            pair.atom0 = aiInner.nextAtom();
            if (needUpdate) {
        		needUpdate = false;
    			pair.atom1 = targetAtom;
            }
        }
        if(!aiInner.hasNext()) {
        	advanceLists();
        	needUpdate = true;
        }
        return pair;
    }
    
    public AtomSet peek() {
        pair.atom1 = (Atom)aiInner.peek();
        return pair;
    }
    
    public void unset() {
        aiInner.unset();
    }

    /**
     * Returns 2, indicating that this is a pair iterator.
     */
    public int nBody() {
        return 2;
    }
    
    public void reset() {
        if(targetAtom == null) {
            unset();
            return;
        }
        inCentralCell = true;
        upListNow = (direction != IteratorDirective.Direction.DOWN);
        neighborIterator.checkDimensions();
        lattice.latticeIndex(centralCell.latticeArrayIndex,latticeIndex);
        neighborIterator.setSite(latticeIndex);
        neighborIterator.setDirection(upListNow ? IteratorDirective.Direction.UP : IteratorDirective.Direction.DOWN);
        neighborIterator.reset();
        needUpdate = false;
        
        //start with targetMolecule's cell
        atomToArrayListFixed.setArrayList(centralCell.occupants());
        if (upListNow) {
            aiSeqDirectableUp.setAtom(targetAtom);
            aiSeqDirectableUp.reset();
            pair.atom0 = targetAtom;
            if (aiSeqDirectableUp.hasNext()) {
                aiInner = aiSeqDirectableUp;
                return;
            }
        }
        if (doGoDown) {
            aiSeqDirectableDn.setAtom(targetAtom);
            aiSeqDirectableDn.reset();
            pair.atom1 = targetAtom;
            if (aiSeqDirectableDn.hasNext()) {
                upListNow = false;
                aiInner = aiSeqDirectableDn;
                return;
            }
        }
        inCentralCell = false;
        advanceLists();
    }
    
    /**
     * Indicates allowed direction for iteration, relative to specified target
     * atom. Specification of a null direction indicates iteration in both directions
     * relative to the target. Direction is determined by ordering within occupant
     * list of cell of target atom, and then by the cell ordering of neighboring cells.
     */
    public void setDirection(Direction direction) {
        this.direction = direction;
        doGoDown = (direction != IteratorDirective.Direction.UP);
        neighborIterator.setDirection(direction);
        nbrCellListIteratorInner.setDirection(direction);
    }

    /**
     * Sets the target molecule with which all pairs are formed.  Molecule
     * is determined from the first atom of the array, which may be the molecule
     * itself or an atom that is part of it.  If the atom is null or is not 
     * in one of the species given at construction, no iterates will be returned.
     * @throws NullPointerException
     *          if targetAtoms is null; use AtomSet.NULL instead
     * @throws IllegalArgumentException
     *          if targetAtoms.count() is not 0 or 1
     */
    public void setTarget(AtomSet targetAtoms) {
        switch(targetAtoms.count()) {
        case 0: 
            targetAtom = null;
            break;
        case 1:
            targetAtom = targetAtoms.getAtom(0);
            break;
        default:
            throw new IllegalArgumentException("Can specify at most one target atom to iterator");
        }
    }
    
    public void setCentralCell(Cell cell) {
        centralCell = cell;
    }

    // Moves to next neighbor-cell list that can provide an iterate
    // This should be invoked only if aiInner.hasNext is false
    private void advanceLists() {
        if (inCentralCell && upListNow && doGoDown) {
            aiSeqDirectableDn.setAtom(targetAtom);
            aiSeqDirectableDn.reset();
            if (aiSeqDirectableDn.hasNext()) {
                upListNow = false;
                aiInner = aiSeqDirectableDn;
                return;
            }
        }
        if (direction == null && inCentralCell) {
            upListNow = true;
        }
        inCentralCell = false;
        aiInner = aiSeq;//need to switch from aiSeqDirectableXX
        do {
            //advance neighbor cell 
            if(neighborIterator.hasNext()) {
                aiSeq.setList(((Cell)neighborIterator.next()).occupants());
                aiSeq.reset();
            } else if (upListNow && doGoDown) {
                upListNow = false;
                neighborIterator.setDirection(IteratorDirective.Direction.DOWN);
                neighborIterator.reset();
                aiSeq.setList(((Cell)neighborIterator.next()).occupants());
                aiSeq.reset();
            } else {//no more cells
                break;
            }
        } while(!aiInner.hasNext());
    }//end of advanceCell

    /**
     * @return Returns the cellIterator.
     */
    public CellLattice.NeighborIterator getNbrCellIterator() {
        return neighborIterator;
    }
   
    private final ApiSequence1A nbrCellListIteratorInner;//used only by allAtoms
    private final ApiInnerFixed nbrCellListIteratorUp;//used only by allAtoms
    private final ApiInnerFixed nbrCellListIteratorDn;//used only by allAtoms
    private final CellLattice.NeighborIterator neighborIterator;
    private final AtomIteratorArrayList aiSeqDirectableUp, aiSeqDirectableDn;
    private final AtomIteratorArrayListSimple aiSeq;
    private final AtomToArrayListFixed atomToArrayListFixed;
    private final AtomIteratorSinglet aiOuter;
    private final AtomPair pair = new AtomPair();
    private final int[] latticeIndex;
    private IteratorDirective.Direction direction;
    private boolean doGoDown, upListNow;
    private boolean inCentralCell;
    private boolean needUpdate;
    private Atom targetAtom;
    private Cell centralCell;
    private final PhaseAgentManager phaseAgentManager;
    
    private CellLattice lattice;
    
    private AtomIterator aiInner;

}
