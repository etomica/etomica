/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.ApiInnerFixed;
import etomica.Atom;
import etomica.AtomIteratorListSimple;
import etomica.AtomIteratorSinglet;
import etomica.AtomList;
import etomica.AtomTreeNode;
import etomica.AtomTreeNodeGroup;
import etomica.AtomsetIteratorMolecule;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Species;
import etomica.IteratorDirective.Direction;
import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.action.AtomsetDetect;
import etomica.lattice.CellLattice;

/**
 * Gives pairs formed from the molecules of two different species in a phase,
 * taking one molecule of one species with all neighboring molecules of the other.
 * Species are specified at construction and cannot be changed afterwards.  The
 * 1-molecule species is identified via the setTarget method, and may be changed
 * from one use of the iterator to the next.  Iteration is performed using cell lists,
 * which defines the neighboring molecules.  Direction is related to ordering of
 * species, and (unlike ApiIntraspecies1ACell) is not connected to cell ordering.
 */

public class ApiInterspecies1ACell implements AtomsetIteratorMolecule, AtomsetIteratorCellular {

	/**
	 * Constructor makes iterator that must have phase specified and then be 
	 * reset() before iteration.
     * 
     * @param D the dimension of the space of the simulation (used to construct cell iterators)
     * @param species length = 2 array with the (different) species whose molecules are interacting 
     */
	public ApiInterspecies1ACell(int D, Species[] species) {
        if(species[0] == null || species[1] == null) throw new NullPointerException("Constructor of ApiInterspecies1ACell requires two non-null species");
        if(species[0] == species[1]) throw new IllegalArgumentException("Constructor of ApiInterspecies1ACell requires two different species");
        //arrange so that species0 preceeds species1
        if(species[0].getIndex() < species[1].getIndex()) {
            species0 = species[0];
            species1 = species[1];
        } else {
            species0 = species[1];
            species1 = species[0];
        }
        index0 = species0.getIndex();
        index1 = species1.getIndex();

        neighborIterator = new CellLattice.NeighborIterator(D);
        aiOuter = new AtomIteratorSinglet();
        aiInner = new AtomIteratorListSimple();
        listIterator = new ApiInnerFixed(aiOuter, aiInner);//used only by allAtoms

        // cell iterator always iterates both up and down the cells, 
        // regardless of what happens with this.setDirection, because the
        // iteration direction for different species is determined by the
        // species order; setDirection for this iterator is considered
        // in connection to the target and the species order.
        // in contrast, for intraspecies iteration direction is related to the cell ordering
        neighborIterator.setDirection(null);
        setPhase(null);
	}

	public void setPhase(Phase phase) {
        this.phase = phase;
        if(phase != null) {
            neighborIterator.setLattice(phase.getLattice());
            agentNode0 = (AtomTreeNodeGroup)phase.getAgent(species0).node;
            agentNode1 = (AtomTreeNodeGroup)phase.getAgent(species1).node;
        }
        identifyTargetMolecule();
	}
    
    /**
     * Performs action on all iterates.
     */
    public void allAtoms(AtomsetAction action) {
        if(pair[0] == null) return;
        aiOuter.setAtom(pair[0]);
        neighborIterator.checkDimensions();
        NeighborCell cell = ((AtomSequencerCell)targetMolecule.seq).cell;
        int[] index = lattice.latticeIndex(cell.latticeArrayIndex);
        
        //get pairs in targetMolecule's cell
        AtomList list = cell.occupants()[innerIndex];
        aiInner.setList(list);
        listIterator.allAtoms(action);

        //loop over neighbor cells
        neighborIterator.setSite(index);
        neighborIterator.reset();
        while(neighborIterator.hasNext()) {
            NeighborCell neighborCell = (NeighborCell)neighborIterator.next(); 
            aiInner.setList(neighborCell.occupants()[innerIndex]);
            if(aiInner.size() > 0) listIterator.allAtoms(action);
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
	public boolean contains(Atom[] atoms) {
        if(atoms==null || atoms[0]==null || atoms[1]==null || atoms[0]==atoms[1]) return false;
        AtomsetDetect detector = new AtomsetDetect(atoms);
        allAtoms(detector);
        return detector.detectedAtom();
	}

    public boolean hasNext() {
        return aiInner.hasNext();
    }
    
    public Atom[] next() {
        if(!hasNext()) return null;
        pair[1] = aiInner.nextAtom();
        if(!aiInner.hasNext()) {
            advanceLists();
        }
        return pair;
    }
    
    public Atom[] peek() {
        pair[1] = aiInner.peek()[0];
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
        if(pair[0] == null) {
            unset();
            return;
        }
        neighborIterator.checkDimensions();
        NeighborCell cell = ((AtomSequencerCell)targetMolecule.seq).cell;
        int[] index = lattice.latticeIndex(cell.latticeArrayIndex);
        neighborIterator.setSite(index);
        neighborIterator.reset();
        
        //start with targetMolecule's cell
        aiInner.setList(cell.occupants()[innerIndex]);
        aiInner.reset();

        if(!aiInner.hasNext()) { 
            advanceLists();
        }
    }
    
    /**
     * Indicates allowed direction for iteration, relative to specified target
     * atom. If the specified direction is consisent with the direction from the
     * target species to the non-target species (as given by their species index --
     * UP is direction from smaller index to larger index) direction, iteration
     * is performed; if specified direction contradicts species direction, no
     * iteration is performed.  Specification of a null direction indicates no
     * limitation, and iteration will be performed if a legitimate target atom 
     * is specified. 
     */
    public void setDirection(Direction direction) {
        this.direction = direction;
        setupIterators();
    }

    /**
     * Sets the target molecule with which all pairs are formed.  Molecule
     * is determined from the first atom of the array, which may be the molecule
     * itself or an atom that is part of it.  If the atom is null or is not 
     * in one of the species given at construction, no iterates will be returned.
     */
    public void setTarget(Atom[] targetAtoms) {
        switch(targetAtoms.length) {
            case 0: 
                targetAtom = null;
                break;
            case 1: 
                targetAtom = targetAtoms[0];
                break;
            default:
                if(targetAtoms[1] != null) throw new IllegalArgumentException("Specification of more than two target atoms is not supported");
        }
        identifyTargetMolecule();
    }

    
    // Moves to next neighbor-cell list that can provide an iterate
    // This should be invoked only if aiInner.hasNext is false
    private void advanceLists() {
        do {
              //advance neighbor cell 
            if(neighborIterator.hasNext()) {
                aiInner.setList(((NeighborCell)neighborIterator.next()).occupants()[innerIndex]);
                aiInner.reset();
            } else {//no more cells
                break;
            }
        } while(!aiInner.hasNext());
    }//end of advanceCell
    
    /**
     * Finds target molecule as indicated by the target atom.  Sets
     * target molecule to null if target atom is null, phase is null, or
     * atom is not part of either species.
     */
    private void identifyTargetMolecule() {
        if(targetAtom == null || phase == null) {
            targetMolecule = null;
        } else {
            AtomTreeNode targetNode = targetAtom.node.childWhereDescendedFrom(agentNode0);

            if(targetNode != null) {    //target is species0, loop over species1
                allowedDirection = IteratorDirective.UP;
                targetMolecule = targetNode.atom();
                innerIndex = index1;
                
            } else {                    //target is not species0
                targetNode = targetAtom.node.childWhereDescendedFrom(agentNode1);
                
                if(targetNode != null) {//target is species1, loop over species0
                    allowedDirection = IteratorDirective.DOWN;
                    targetMolecule = targetNode.atom();
                    innerIndex = index0;
                    
                } else {                //target not in either species
                    targetMolecule = null;
                }
            }
        }
        setupIterators();
    }
    
    /**
     * Completes setup of iterators, checking that specified direction
     * is consistent with target and species ordering. Leaves this
     * iterator unset.
     */
    private void setupIterators() {
        if(direction == null || direction == allowedDirection) {
            pair[0] = targetMolecule;
            aiOuter.setAtom(targetMolecule);//targetMolecule may be null here
        } else {
            pair[0] = null;
            aiOuter.setAtom(null);
        }
        unset();
    }

    /**
     * @return Returns the cellIterator.
     */
    public CellLattice.NeighborIterator getNbrCellIterator() {
        return neighborIterator;
    }
   

    private final ApiInnerFixed listIterator;//used only by allAtoms
    private final CellLattice.NeighborIterator neighborIterator;
    private final int index0, index1;
    private final AtomIteratorListSimple aiInner;
    private final AtomIteratorSinglet aiOuter;
    
    private final Atom[] pair = new Atom[2];
    private int innerIndex;
    
    private final Species species0, species1;
    private AtomTreeNodeGroup agentNode0, agentNode1;
    private Atom targetMolecule, targetAtom;
    private Phase phase;
    private IteratorDirective.Direction allowedDirection, direction;
    private CellLattice lattice;
}
