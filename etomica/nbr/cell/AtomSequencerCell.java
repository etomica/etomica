/*
 * History
 * Created on Nov 23, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.Atom;
import etomica.AtomIteratorListSimple;
import etomica.AtomLinker;
import etomica.AtomReservoir;
import etomica.AtomSequencer;
import etomica.AtomTreeNodeGroup;
import etomica.nbr.cell.IteratorFactoryCell.CellSequencer;


/**
 * Sequencer used for atoms being cell listed.
 */
public final class AtomSequencerCell extends AtomSequencer implements CellSequencer {
    
    public NeighborCell cell;         //cell currently occupied by this atom
    public NeighborCellManager neighborCellManager;    //cell lattice in the phase occupied by this atom
    private int listIndex;
    public final AtomLinker nbrLink;  //linker used to arrange atom in sequence according to cells
    
    public AtomSequencerCell(Atom a) {
        super(a);
        nbrLink = new AtomLinker(a);
    }
    
    //CellSequencer interface method
    public void latticeChangeNotify() {
        this.assignCell();
        if(atom.node.isLeaf()) return;
        else {
            AtomIteratorListSimple iterator = new AtomIteratorListSimple(((AtomTreeNodeGroup)atom.node).childList);
            while(iterator.hasNext()) ((CellSequencer)atom.seq).latticeChangeNotify();
        }
    }

    //CellSequencer interface method
    public NeighborCell cell() {return cell;}

    public int listIndex() {return listIndex;}

    public void remove() {
        super.remove();
        nbrLink.remove();
    }
        
    public void addBefore(AtomLinker newNext) {
        //newNext will never be one of the cell tabs
        super.addBefore(newNext);
        if(neighborCellManager != null) assignCell();
        else nbrLink.remove();
        
    /*    while(newNext.atom == null && newNext != this) newNext = newNext.next;
        if(newNext == this) {//this is the first and only non-tab entry in the list
            nextFixed = previousFixed = this;
            return;
        }
        nextFixed = (Sequencer)newNext;
        previousFixed = nextFixed.previousFixed;
        previousFixed.nextFixed = this;
        nextFixed.previousFixed = this;*/
	}
	/**
	 * Reshuffles position of "neighbor" links without altering the regular links.
	 */
	public void moveBefore(AtomLinker newNext) {
	    nbrLink.moveBefore(newNext);
	}


    /**
     * Returns true if this atom preceeds the given atom in the atom sequence.
     * Returns false if the given atom is this atom, or (of course) if the
     * given atom instead preceeds this one.
     */
     //this method needs to be fixed
    public boolean preceeds(Atom a) {
        throw new RuntimeException("IteratorFactoryCell.CellSequencer.preceeds method not yet implemented");
        //want to return false if atoms are the same atoms
      /*  if(a == null) return true;
        if(atom.node.parentGroup() == a.node.parentGroup()) {
            if(((Sequencer)atom.seq).site().equals(cell)) {
                //this isn't correct
                return atom.node.index() < a.node.index();//works also if both parentGroups are null
            }
            else return ((Sequencer)atom.seq).site().preceeds(cell);
        }
        int thisDepth = atom.node.depth();
        int atomDepth = a.node.depth();
        if(thisDepth == atomDepth) return atom.node.parentGroup().seq.preceeds(a.node.parentGroup());
        else if(thisDepth < atomDepth) return this.preceeds(a.node.parentGroup());
        else /*if(this.depth > atom.depth)* / return atom.node.parentGroup().seq.preceeds(a);
        */
    }
    
    /**
     * Method called when the parent of the atom is changed.
     * By the time this method is called, the atom has been placed
     * in the childList of the given parent (if it is not null).
     */
    public void setParentNotify(AtomTreeNodeGroup newParent) {
        if(newParent == null || newParent instanceof AtomReservoir.ReservoirAtomTreeNode) {
            cell = null;
            neighborCellManager = null;
            nbrLink.remove();//08/12/03 (DAK) added this line
            return;
        }
        //get cell lattice for the phase containing the parent
        neighborCellManager = newParent.parentPhase().getLattice();
        listIndex = newParent.parentSpeciesAgent().node.index();
        while (listIndex+1 > neighborCellManager.getListCount()) {
            neighborCellManager.addList();
        }
        assignCell();
    }

//Determines appropriate cell and assigns it
    public void assignCell() {
        int[] latticeIndex = neighborCellManager.getCellLattice().getPrimitive().latticeIndex(atom.coord.position(), neighborCellManager.getCellLattice().getDimensions());
        NeighborCell newCell = (NeighborCell)neighborCellManager.getCellLattice().site(latticeIndex);
        if(newCell != cell) {assignCell(newCell);}
    }
//Assigns atom to given cell
    public void assignCell(NeighborCell newCell) {
        cell = newCell;
        if(newCell == null) {
        	nbrLink.remove();//08/12/03 (DAK) added this line
        } else {
	        this.moveBefore(((AtomLinker.Tab[])newCell.agents[0])[listIndex].nextTab);
        }
    }//end of assignCell
    
    public static final AtomSequencer.Factory FACTORY = new AtomSequencer.Factory() {
        public AtomSequencer makeSequencer(Atom atom) {return new AtomSequencerCell(atom);}
        public Class sequencerClass() {return AtomSequencerCell.class;}
    };
}