/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.action.AtomsetDetect;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.iterator.ApiInterArrayList;
import etomica.atom.iterator.ApiIntraArrayList;
import etomica.atom.iterator.AtomPairIterator;
import etomica.atom.iterator.IteratorDirective;
import etomica.lattice.CellLattice;
import etomica.lattice.RectangularLattice;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;
import etomica.space.BoundaryPeriodic;

/**
 * Returns iterates formed from all cell-based neighbor pairs.
 */

public class ApiAACell implements AtomPairIterator, AtomsetIteratorCellular, java.io.Serializable {

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
     */
	public ApiAACell(int D, double range, PhaseAgentManager agentManager) {
        cellIterator = new RectangularLattice.Iterator(D);
        neighborIterator = new CellLattice.NeighborIterator(D, range);
        neighborIterator.setDirection(IteratorDirective.Direction.UP);
        interListIterator = new ApiInterArrayList();
        intraListIterator = new ApiIntraArrayList();
        listIterator = intraListIterator;
        phaseAgentManager = agentManager;
	}

	public void setPhase(Phase phase) {
        if(this.phase == phase) {
            unset();
            return;
        }
        this.phase = phase;
        NeighborCellManager[] cellManagers = (NeighborCellManager[])phaseAgentManager.getAgents();
        CellLattice lattice = cellManagers[phase.getIndex()].getLattice();
        cellIterator.setLattice(lattice);
		neighborIterator.setLattice(lattice);
        neighborIterator.setPeriod(phase.getBoundary().getDimensions());
        neighborIterator.setPeriodicity(((BoundaryPeriodic)phase.getBoundary()).getPeriodicity());
        unset();
	}
    
    /**
     * Performs action on all iterates.
     */
    public void allAtoms(AtomsetAction action) {
        cellIterator.reset();
        neighborIterator.checkDimensions();
        while(cellIterator.hasNext()) {//outer loop over all cells
            //get cell without advancing -- advance is done via nextIndex,
            // below
            Cell cell = (Cell)cellIterator.peek();
            AtomArrayList list = cell.occupants();
            
            //no molecules of species in cell
            if(list.isEmpty()) {
                cellIterator.nextIndex();
                continue;
            }
            
            //consider pairs formed from molecules in cell
            intraListIterator.setList(list);
            intraListIterator.allAtoms(action);

            //loop over neighbor cells
            interListIterator.setOuterList(list);
            neighborIterator.setSite(cellIterator.nextIndex());
            neighborIterator.reset();
            while(neighborIterator.hasNext()) {
                Cell neighborCell = (Cell)neighborIterator.next(); 
                interListIterator.setInnerList(neighborCell.occupants());
                if(neighborCell.occupants().size() > 0) interListIterator.allAtoms(action);
            }
        }//end of outer loop over cells
    }//end of allAtoms
    
	/**
     * Returns the number of atom pairs the iterator will return if reset and
     * iterated in its present state.
     */
	public int size() {
        AtomsetCount counter = new AtomsetCount();
        allAtoms(counter);
        return counter.callCount();
	}
	
	/**
     * Indicates whether the given atom pair will be among the iterates given by
     * the iterator if reset in its present state. True only if an iterated pair
     * would match the atoms as ordered in the given array.
     */
	public boolean contains(AtomSet atoms) {
        if(!(atoms instanceof AtomPair) || ((AtomPair)atoms).atom0 == ((AtomPair)atoms).atom1) return false;
        AtomsetDetect detector = new AtomsetDetect(atoms);
        allAtoms(detector);
        return detector.detectedAtom();
	}

    public boolean hasNext() {
        return listIterator.hasNext();
    }
    
    public AtomSet next() {
        return nextPair();
    }
    
    public AtomPair nextPair() {
        if(!hasNext()) return null;
        AtomPair nextPair = listIterator.nextPair();
        nextPair.copyTo(pair);
        if(!listIterator.hasNext()) {
            advanceLists();
        }
        return pair;
    }
    
    public AtomSet peek() {
        return listIterator.peek();
    }
    
    public void unset() {
        listIterator.unset();
    }

    /**
     * Returns 2, indicating that this is a pair iterator.
     */
    public int nBody() {
        return 2;
    }
    
    public void reset() {
        cellIterator.reset();
        neighborIterator.checkDimensions();
        neighborIterator.unset();
        listIterator.unset();
        advanceLists();
        //System.out.println("reset in ApiAACell");
    }//end of reset
    
    // Moves to next pair of lists that can provide an iterate
    // This should be invoked only if listIterator.hasNext is false
    private void advanceLists() {
        do {
              //advance neighbor cell
            if(neighborIterator.hasNext()) {
                interListIterator.setInnerList(((Cell)neighborIterator.next()).occupants());
                listIterator = interListIterator;
                interListIterator.reset();

                //advance central cell and set up neighbor cell iterator if
                // central cell has some molecules
            } else if(cellIterator.hasNext()) {
                AtomArrayList list = ((Cell)cellIterator.peek()).occupants();
                neighborIterator.setSite(cellIterator.nextIndex());

                if(!list.isEmpty()) {//central cell has molecules
                    interListIterator.setOuterList(list); //for neighbor-cell looping
                    intraListIterator.setList(list);//for intra-cell looping
                    neighborIterator.reset();

                    listIterator = intraListIterator;
                    intraListIterator.reset();
                        
                } else {//no molecules in central cell
                    neighborIterator.unset();
                    listIterator.unset();
                }
            } else {//no more cells at all
                break;
            }
        } while(!listIterator.hasNext());
    }//end of advanceCell
    
    /**
     * @return Returns the cellIterator.
     */
    public CellLattice.NeighborIterator getNbrCellIterator() {
        return neighborIterator;
    }
   
    private AtomPairIterator listIterator;
    private Phase phase;
    private final ApiIntraArrayList intraListIterator;
    private final ApiInterArrayList interListIterator;
    private final CellLattice.NeighborIterator neighborIterator;
    private final RectangularLattice.Iterator cellIterator;
    private final PhaseAgentManager phaseAgentManager;

    private final AtomPair pair = new AtomPair();
}
