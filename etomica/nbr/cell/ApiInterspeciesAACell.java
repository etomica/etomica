/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.AtomPair;
import etomica.AtomPairIterator;
import etomica.AtomSet;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Species;
import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.action.AtomsetDetect;
import etomica.atom.AtomList;
import etomica.atom.AtomPairVector;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.ApiInnerFixed;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.atom.iterator.AtomsetIteratorPhaseDependent;
import etomica.lattice.CellLattice;
import etomica.lattice.RectangularLattice;

/**
 * Returns iterates formed from all molecule pairs of two species. Looping is
 * such that each pair of cells is considered once. For each pair (cellA,
 * cellB), iteration is performed using the species0 molecules in cellA with
 * species1 in cellB, and species1 in cellA with species0 in cellB.
 */

public class ApiInterspeciesAACell implements AtomsetIteratorPhaseDependent, AtomsetIteratorCellular, 
        AtomPairIterator {

	/**
	 * Constructor makes iterator that must have phase specified and then be 
	 * reset() before iteration.
     * 
     * @param D the dimension of the space of the simulation (used to construct cell iterators)
     * @param species length = 2 array with the (different) species whose molecules are interacting 
     */
	public ApiInterspeciesAACell(int D, Species[] species) {
        cellIterator = new RectangularLattice.Iterator(D);
        neighborIterator = new CellLattice.NeighborIterator(D);
        neighborIterator.setDirection(IteratorDirective.UP);
        listIterator = ApiBuilder.makeInterlistIterator();
        aiInner = ((AtomIteratorListSimple)listIterator.getInnerIterator());
        aiOuter = ((AtomIteratorListSimple)listIterator.getOuterIterator());
        index0 = species[0].getIndex()-1;
        index1 = species[1].getIndex()-1;
        if(index0 == index1) throw new IllegalArgumentException("Intergroup iterator cannot be constructed with a single species");
	}

	public void setPhase(Phase phase) {
        cellIterator.setLattice(phase.getLattice());
		neighborIterator.setLattice(phase.getLattice());
        neighborIterator.setPeriod(phase.boundary().dimensions());
        unset();
	}

    /**
     * Performs action on all iterates.
     */
    public void allAtoms(AtomsetAction action) {
        cellIterator.reset();
        neighborIterator.checkDimensions();
        while(cellIterator.hasNext()) {//outer loop over all cells
            //get cell without advancing -- advance is done via nextIndex, below
            NeighborCell cell = (NeighborCell)cellIterator.peek();
            AtomList list0 = cell.occupants()[index0];
            AtomList list1 = cell.occupants()[index1];
            boolean do0 = !list0.isEmpty();
            boolean do1 = !list1.isEmpty();
            
            //no molecules of either species in cell
            if(!do0 && !do1) {
                cellIterator.nextIndex();
                continue;
            }
            
            //consider pairs formed from molecules in cell
            if(do0 && do1) {
                aiOuter.setList(list0);
                aiInner.setList(list1);
                listIterator.allAtoms(action);
            }

            //loop over neighbor cells
            neighborIterator.setSite(cellIterator.nextIndex());
            if(do0) neighborCellLoop(list0, index1, action);
            if(do1) neighborCellLoop(list1, index0, action);
            
            //alternative structure that loops over neighbor cells only once
//            neighborIterator.reset();
//            while(neighborIterator.hasNext()) {//inner loop over neighbor cells
//                NeighborCell neighborCell = (NeighborCell)neighborIterator.next(); 
//                if(do0) {
//                    aiOuter.setList(list0);
//                    aiInner.setList(neighborCell.occupants()[index1]);
//                    if(aiInner.size() > 0) listIterator.allAtoms(action);
//                }
//                if(do1) {
//                    aiOuter.setList(list1);
//                    aiInner.setList(neighborCell.occupants()[index0]);
//                    if(aiInner.size() > 0) listIterator.allAtoms(action);
//                }
//            }//end of inner loop over neighbor cells
            
            
        }//end of outer loop over cells
    }//end of allAtoms
    
    //used by allAtoms
    private void neighborCellLoop(AtomList list, int index, AtomsetAction action) {
        aiOuter.setList(list);
        neighborIterator.reset();
        while(neighborIterator.hasNext()) {
            NeighborCell neighborCell = (NeighborCell)neighborIterator.next(); 
            aiInner.setList(neighborCell.occupants()[index]);
            if(aiInner.size() > 0) listIterator.allAtoms(action);
        }
    }
	
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
        return listIterator.hasNext();
    }
    
    public AtomSet next() {
        return nextPair();
    }
    
    public AtomPair nextPair() {
        if(!hasNext()) return null;
        AtomPair nextPair = listIterator.nextPair();//returns AtomPair
        nextPair.copyTo(pair);//copy to AtomPairVector
        pair.nearestImageVector = neighborIterator.getNearestImageVector();
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
        doListSwap = false;
        listIterator.unset();
        advanceLists();
    }
    
    // Moves to next pair of lists that can provide an iterate
    // This should be invoked only if listIterator.hasNext is false
    private void advanceLists() {
        do {
              //advance neighbor cell 
            if(neighborIterator.hasNext()) {
                aiInner.setList(((NeighborCell)neighborIterator.next()).occupants()[innerIndex]);
                listIterator.reset();

             //change from cellA/index0-cellB/index1 iteration to cellA/index1-cellB/index0
            } else if(doListSwap) {
                //we'll get here only if list1 is not empty
                aiOuter.setList(centralCell.occupants()[index1]);
                innerIndex = index0;
                neighborIterator.reset();
                doListSwap = false;

                //advance central cell, prepare to iterate over pairs in it, and set up neighbor cell iterator
            } else if(cellIterator.hasNext()) {
                centralCell = (NeighborCell)cellIterator.peek();
                neighborIterator.setSite(cellIterator.nextIndex());
 
                AtomList list0 = centralCell.occupants()[index0];
                AtomList list1 = centralCell.occupants()[index1];

                if(!list0.isEmpty()) {//central cell has (at least) molecules of species0 
                    aiOuter.setList(list0);
                    innerIndex = index1;
                    neighborIterator.reset();

                    if(!list1.isEmpty()) {//central cell has molecules of both species
                        aiInner.setList(list1);
                        listIterator.reset();//this is only place where listIterator is set with both centralCell lists
                        doListSwap = true;//swap because want centralCell to be outer loop for both species
                        
                    } else {//central cell has molecules of species0 only
                        listIterator.unset();
                        doListSwap = false;
                    }
                    
                } else if(!list1.isEmpty()) {//central cell has molecules of species1 only
                    aiOuter.setList(list1);
                    innerIndex = index0;
                    neighborIterator.reset();
                    listIterator.unset();
                    doListSwap = false;
                    
                } else {//no molecules of either species in central cell
                    neighborIterator.unset();
                    doListSwap = false;
                    listIterator.unset();
                    advanceLists();
                }
            } else {//no more cells at all
                break;
            }
            
            //alternative approach that does a single iteration of neighbor cells
            
            //consider pairs formed from cellA/index1 with cellB/index0
//            if(doListSwap) {
//                aiOuter.setList(centralCell.occupants()[index1]);
//                aiInner.setList(neighborCell.occupants()[index0]);
//                listIterator.reset();
//                doListSwap = false;
//                
//              //advance neighbor cell and consider pairs formed from cellA/index0 with cellB/index1
//            } else if(neighborIterator.hasNext()) {
//                neighborCell = (NeighborCell)neighborIterator.next();
//                aiOuter.setList(centralCell.occupants()[index0]);
//                aiInner.setList(neighborCell.occupants()[index1]);
//                listIterator.reset();
//                doListSwap = true;
//            }
//            //advance central cell and set up neighbor cell iterator if central cell has some molecules
//        } else if(cellIterator.hasNext()) {
//            centralCell = (NeighborCell)cellIterator.peek();
//            neighborIterator.setSite(cellIterator.nextIndex());
//            AtomList[] occupantLists = centralCell.occupants();
//            //do nbrCell iteration only if central cell has occupants
//            if(occupantLists[index0].size() > 0 || occupantLists[index1].size() > 0) {
//                neighborIterator.reset();
//            }
//        } else {//no more cells at all
//            break;
//        }
                

        } while(!listIterator.hasNext());
    }//end of advanceCell
    
    /**
     * @return Returns the cellIterator.
     */
    public CellLattice.NeighborIterator getNbrCellIterator() {
        return neighborIterator;
    }
   
    private final ApiInnerFixed listIterator;
    private final CellLattice.NeighborIterator neighborIterator;
    private final RectangularLattice.Iterator cellIterator;
    private final int index0, index1;
    private final AtomIteratorListSimple aiInner, aiOuter;
    
    private NeighborCell centralCell, neighborCell;
    private boolean doListSwap;
    private final AtomPairVector pair = new AtomPairVector();
    private int innerIndex;
}
