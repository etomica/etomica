/*
 * History
 * Created on Nov 23, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.Atom;
import etomica.AtomLinker;
import etomica.AtomSequencerFactory;


/**
 * Sequencer used for atoms being cell listed.
 */
public class AtomSequencerCell extends AtomLinker implements Cellular {
    
    NeighborCell cell;       //cell currently occupied by this atom
    final AtomLinker nbrLink;  //linker used to arrange atom in sequence within the cells
    
    public AtomSequencerCell(Atom a) {
        super(a);
        nbrLink = new AtomLinker(a);
    }
    
    //Cellular interface method
    public NeighborCell getCell() {return cell;}

//    public int listIndex() {return listIndex;}

    public void remove() {
        super.remove();
        nbrLink.remove();
    }
    
//    /**
//     * Method called when the parent of the atom is changed.
//     * By the time this method is called, the atom has been placed
//     * in the childList of the given parent (if it is not null).
//     */
//    public void setParentNotify(AtomTreeNodeGroup newParent) {
//        if(newParent == null || newParent instanceof AtomReservoir.ReservoirAtomTreeNode) {
//            cell = null;
//            neighborCellManager = null;
//            nbrLink.remove();//08/12/03 (DAK) added this line
//            return;
//        }
//        //get cell lattice for the phase containing the parent
//        neighborCellManager = newParent.parentPhase().getLattice();
//        listIndex = newParent.parentSpeciesAgent().node.index();
//        while (listIndex+1 > neighborCellManager.getListCount()) {
//            neighborCellManager.addList();
//        }
//        assignCell();
//    }

    public static final AtomSequencerFactory FACTORY = new AtomSequencerFactory() {
        public AtomLinker makeSequencer(Atom atom) {return new AtomSequencerCell(atom);}
    };
}