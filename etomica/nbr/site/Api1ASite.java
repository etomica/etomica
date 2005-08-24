/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica.nbr.site;

import etomica.Phase;
import etomica.IteratorDirective.Direction;
import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.action.AtomsetDetect;
import etomica.atom.AtomPair;
import etomica.atom.AtomPairIterator;
import etomica.atom.AtomSet;
import etomica.atom.iterator.AtomsetIteratorMolecule;
import etomica.lattice.CellLattice;
import etomica.lattice.RectangularLatticeNbrIterator;
import etomica.lattice.RectangularLatticeNbrIteratorAdjacent;
import etomica.nbr.cell.AtomSequencerCell;
import etomica.nbr.cell.Cell;
import etomica.space.BoundaryPeriodic;

/**
 * Gives pairs formed from the molecules of a species in a phase, taking one
 * molecule the species with all of its other neighboring molecules. Species is
 * specified at construction and cannot be changed afterwards. Iteration is
 * performed using cell lists, which defines the neighboring molecules.
 * Direction is related to ordering of cells and, within a cell, ordering of
 * molecules in cell's occupant list.
 */

public class Api1ASite implements AtomsetIteratorMolecule, AtomPairIterator {
    
	/**
	 * Constructor makes iterator that must have phase specified and then be 
	 * reset() before iteration.
     * 
     * @param D the dimension of the space of the simulation (used to construct cell iterators)
     * @param species length = 2 array with the (different) species whose molecules are interacting 
     */
	public Api1ASite(int D) {
        neighborIterator = new RectangularLatticeNbrIteratorAdjacent(D);
        
        latticeIndex = new int[D];

        neighborIterator.setDirection(null);
        setPhase(null);
	}

	public void setPhase(Phase phase) {
        if(phase != null) {
            lattice = (CellLattice)phase.getLattice();
            neighborIterator.setLattice(lattice);
            neighborIterator.setPeriodicity(((BoundaryPeriodic)phase.boundary()).getPeriodicity());
        }
	}

    /**
     * Performs action on all iterates.
     */
    public void allAtoms(AtomsetAction action) {
        if(pair.atom0 == null) return;
        Cell cell = ((AtomSequencerCell)pair.atom0.seq).getCell();
        lattice.latticeIndex(cell.getLatticeArrayIndex(),latticeIndex);
        
        //loop over neighbor cells
        neighborIterator.setSite(latticeIndex);
        neighborIterator.reset();
        while(neighborIterator.hasNext()) {
            pair.atom1 = ((AtomSite)neighborIterator.next()).getAtom();
            action.actionPerformed(pair);
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
        return neighborIterator.hasNext();
    }
    
    public AtomSet next() {
        return nextPair();
    }
    
    public AtomPair nextPair() {
        if(!hasNext()) return null;
        pair.atom1 = ((AtomSite)neighborIterator.next()).getAtom();
        return pair;
    }
    
    public AtomSet peek() {
        pair.atom1 = ((AtomSite)neighborIterator.peek()).getAtom();
        return pair;
    }
    
    public void unset() {
        neighborIterator.unset();
    }

    /**
     * Returns 2, indicating that this is a pair iterator.
     */
    public int nBody() {
        return 2;
    }
    
    public void reset() {
        if(pair.atom0 == null) {
            unset();
            return;
        }
        AtomSite site = ((AtomSequencerSite)pair.atom0.seq).getSite();
        lattice.latticeIndex(site.latticeArrayIndex,latticeIndex);
        neighborIterator.setSite(latticeIndex);
        neighborIterator.reset();
    }
    
    /**
     * Indicates allowed direction for iteration, relative to specified target
     * atom. Specification of a null direction indicates iteration in both directions
     * relative to the target. Direction is determined by ordering within occupant
     * list of cell of target atom, and then by the cell ordering of neighboring cells.
     */
    public void setDirection(Direction direction) {
        neighborIterator.setDirection(direction);
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
            pair.atom0 = null;
            break;
        case 1:
            pair.atom0 = targetAtoms.getAtom(0);
            break;
        default:
            throw new IllegalArgumentException("Can specify at most one target atom to iterator");
        }
    }

    
    private final RectangularLatticeNbrIterator neighborIterator;
    private final AtomPair pair = new AtomPair();
    private final int[] latticeIndex;
    
    private CellLattice lattice;
    
}
