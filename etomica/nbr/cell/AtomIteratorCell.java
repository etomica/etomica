package etomica.nbr.cell;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.lattice.CellLattice;
import etomica.lattice.RectangularLattice;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;

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
	public AtomIteratorCell(int D, PhaseAgentManager agentManager) {
        cellIterator = new RectangularLattice.Iterator(D);
        atomIterator = new AtomIteratorArrayListSimple();
        phaseAgentManager = agentManager;
	}

	public void setPhase(Phase phase) {
        CellLattice lattice = ((NeighborCellManager)phaseAgentManager.getAgent(phase)).getLattice();
        cellIterator.setLattice(lattice);
        unset();
	}
    
    /**
     * Performs action on all iterates.
     */
    public void allAtoms(AtomsetAction action) {
        cellIterator.reset();
        while(cellIterator.hasNext()) {//outer loop over all cells
            Cell cell = (Cell)cellIterator.next();
            AtomArrayList list = cell.occupants();
            
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
	
    public boolean hasNext() {
        return atomIterator.hasNext();
    }
    
    public final AtomSet next() {
        return nextAtom();
    }
    
    public IAtom nextAtom() {
        if(!hasNext()) return null;
        IAtom nextAtom = atomIterator.nextAtom();
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
                AtomArrayList list = ((Cell)cellIterator.next()).occupants();
                atomIterator.setList(list);
                atomIterator.reset();
            } else {//no more cells at all
                break;
            }
        } while(!atomIterator.hasNext());
    }//end of advanceCell
    
   
    private static final long serialVersionUID = 1L;
    private final AtomIteratorArrayListSimple atomIterator;
    private final RectangularLattice.Iterator cellIterator;
    private final PhaseAgentManager phaseAgentManager;
}
