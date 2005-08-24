/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.Phase;
import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.action.AtomsetDetect;
import etomica.atom.Atom;
import etomica.atom.AtomIterator;
import etomica.atom.AtomList;
import etomica.atom.AtomSet;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.lattice.RectangularLattice;

/**
 * Returns occupants of all cells as iterates.
 */

public class AtomIteratorCell implements AtomIterator, java.io.Serializable {

	/**
     * Constructor makes iterator that must have phase specified and then be
     * reset() before iteration.
     * 
     * @param D
     *            the dimension of the space of the simulation (used to
     *            construct cell iterators)
     */
	public AtomIteratorCell(int D) {
        cellIterator = new RectangularLattice.Iterator(D);
        atomIterator = new AtomIteratorListSimple();
	}

	public void setPhase(Phase phase) {
        cellIterator.setLattice(phase.getLattice());
        unset();
	}
    
    /**
     * Performs action on all iterates.
     */
    public void allAtoms(AtomsetAction action) {
        cellIterator.reset();
        while(cellIterator.hasNext()) {//outer loop over all cells
            Cell cell = (Cell)cellIterator.next();
            AtomList list = cell.occupants();
            
            //consider pairs formed from molecules in cell
            if(!list.isEmpty()) {
                atomIterator.setList(list);
                atomIterator.allAtoms(action);
            }
        }//end of outer loop over cells
    }//end of allAtoms
    
	/**
     * Returns the number of atoms the iterator will return if reset and
     * iterated in its present state.
     */
	public int size() {
        AtomsetCount counter = new AtomsetCount();
        allAtoms(counter);
        return counter.callCount();
	}
	
	/**
     * Indicates whether the given atom will be among the iterates given by
     * the iterator if reset in its present state.
     */
	public boolean contains(AtomSet atoms) {
        if(!(atoms instanceof Atom)) return false;
        AtomsetDetect detector = new AtomsetDetect(atoms);
        allAtoms(detector);
        return detector.detectedAtom();
	}

    public boolean hasNext() {
        return atomIterator.hasNext();
    }
    
    public final AtomSet next() {
        return nextAtom();
    }
    
    public Atom nextAtom() {
        if(!hasNext()) return null;
        Atom nextAtom = atomIterator.nextAtom();
        if(!atomIterator.hasNext()) {
            advanceLists();
        }
        return nextAtom;
    }
    
    public AtomSet peek() {
        return atomIterator.peek();
    }
    
    public void unset() {
        atomIterator.unset();
    }

    /**
     * Returns 1, indicating that this is an atom iterator.
     */
    public int nBody() {
        return 1;
    }
    
    public void reset() {
        cellIterator.reset();
        atomIterator.unset();
        advanceLists();

    }//end of reset
    
    // Moves to next pair of lists that can provide an iterate
    // This should be invoked only if listIterator.hasNext is false
    private void advanceLists() {
        do {
            if(cellIterator.hasNext()) {
                AtomList list = ((Cell)cellIterator.next()).occupants();
                atomIterator.setList(list);
                atomIterator.reset();
            } else {//no more cells at all
                break;
            }
        } while(!atomIterator.hasNext());
    }//end of advanceCell
    
   
    private final AtomIteratorListSimple atomIterator;
    private final RectangularLattice.Iterator cellIterator;

}
