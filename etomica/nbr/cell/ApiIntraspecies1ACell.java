/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.ApiInnerFixed;
import etomica.Atom;
import etomica.AtomIterator;
import etomica.AtomIteratorListSimple;
import etomica.AtomIteratorSequencerList;
import etomica.AtomIteratorSinglet;
import etomica.AtomLinker;
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
 * Gives pairs formed from the molecules of a species in a phase, taking one
 * molecule the species with all of its other neighboring molecules. Species is
 * specified at construction and cannot be changed afterwards. Iteration is
 * performed using cell lists, which defines the neighboring molecules.
 * Direction is related to ordering of cells and, within a cell, ordering of
 * molecules in cell's occupant list.
 */

public class ApiIntraspecies1ACell implements AtomsetIteratorMolecule {

    /**
     * @param species species whose molecules will form the pair iterates
     */
    public ApiIntraspecies1ACell(int D, Species species) {
        this(D, new Species[] {species, species});
    }
    
	/**
	 * Constructor makes iterator that must have phase specified and then be 
	 * reset() before iteration.
     * 
     * @param D the dimension of the space of the simulation (used to construct cell iterators)
     * @param species length = 2 array with the (different) species whose molecules are interacting 
     */
	public ApiIntraspecies1ACell(int D, Species[] species) {
        if(species == null || species.length < 1 || species[0] == null) throw new NullPointerException("Constructor of ApiIntraspecies1A requires two non-null species references to the same instance");
        if(species[0] != species[1]) throw new IllegalArgumentException("Constructor of ApiIntraspecies1A requires references to the same species instance");
        this.species = species[0];
        innerIndex = this.species.getIndex();

        neighborIterator = new CellLattice.NeighborIterator(D);
        aiOuter = new AtomIteratorSinglet();
        aiInnerList = new AtomIteratorListSimple();
        
        //this iterator is used to loop through list of occupants of atoms's cell;
        //construct with AtomToLinker that gives appropriate linker
        aiInnerSeq = new AtomIteratorSequencerList(new AtomIteratorSequencerList.AtomToLinker() {
            public AtomLinker getLinker(Atom atom) {return ((AtomSequencerCell)atom.seq).nbrLink;}
        });
        nbrCellListIterator = new ApiInnerFixed(aiOuter, aiInnerList);//used only by allAtoms
        centralCellListIterator = new ApiInnerFixed(aiOuter, aiInnerSeq);//used only by allAtoms

        aiInnerSeq.setDirection(null);
        aiInnerSeq.setNumToSkip(1);
        neighborIterator.setDirection(null);
        setPhase(null);
	}

	public void setPhase(Phase phase) {
        this.phase = phase;
        if(phase != null) {
            neighborIterator.setLattice(phase.getLattice());
            agentNode = (AtomTreeNodeGroup)phase.getAgent(species).node;
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
        aiInnerSeq.setAtom(pair[0]);
        nbrCellListIterator.allAtoms(action);

        //loop over neighbor cells
        neighborIterator.setSite(index);
        neighborIterator.reset();
        while(neighborIterator.hasNext()) {
            NeighborCell neighborCell = (NeighborCell)neighborIterator.next(); 
            aiInnerList.setList(neighborCell.occupants()[innerIndex]);
            if(aiInnerList.size() > 0) nbrCellListIterator.allAtoms(action);
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
        neighborIterator.setSite(lattice.latticeIndex(cell.latticeArrayIndex));
        neighborIterator.reset();
        
        //start with targetMolecule's cell
        aiInnerSeq.setAtom(pair[0]);
        aiInnerSeq.reset();
        aiInner = aiInnerSeq;

        if(!aiInnerSeq.hasNext()) { 
            advanceLists();
        }
    }
    
    /**
     * Indicates allowed direction for iteration, relative to specified target
     * atom. Specification of a null direction indicates iteration in both directions
     * relative to the target. Direction is determined by ordering within occupant
     * list of cell of target atom, and then by the cell ordering of neighboring cells.
     */
    public void setDirection(Direction direction) {
        aiInnerSeq.setDirection(direction);
        neighborIterator.setDirection(direction);
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
        aiInner = aiInnerList;//need to switch from aiInnerSeq on first call
        do {
              //advance neighbor cell 
            if(neighborIterator.hasNext()) {
                aiInnerList.setList(((NeighborCell)neighborIterator.next()).occupants()[innerIndex]);
                aiInnerList.reset();
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
            AtomTreeNode targetNode = targetAtom.node.childWhereDescendedFrom(agentNode);
            targetMolecule = (targetNode != null) ? targetNode.atom() : null;
        }
        pair[0] = targetMolecule;
        aiOuter.setAtom(targetMolecule);//targetMolecule may be null here
    }

    private final ApiInnerFixed nbrCellListIterator;//used only by allAtoms
    private final ApiInnerFixed centralCellListIterator;//used only by allAtoms
    private final CellLattice.NeighborIterator neighborIterator;
    private final AtomIteratorListSimple aiInnerList;
    private final AtomIteratorSequencerList aiInnerSeq;
    private final AtomIteratorSinglet aiOuter;
    private int innerIndex;
    private final Atom[] pair = new Atom[2];
    
    private final Species species;
    private AtomTreeNodeGroup agentNode;
    private Atom targetMolecule, targetAtom;
    private Phase phase;
    private IteratorDirective.Direction allowedDirection, direction;
    private CellLattice lattice;
    
    private AtomIterator aiInner;
    private boolean finishedCentralCell;

}
