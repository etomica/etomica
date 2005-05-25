/*
 * History
 * Created on Nov 23, 2004 by kofke
 */
package etomica.nbratom.cell;

import etomica.Atom;
import etomica.atom.AtomLinker;
import etomica.atom.AtomSequencerFactory;


/**
 * Sequencer used for atoms being cell listed.  Contains 
 * another linker that is used in the atom lists maintained
 * by the cells.
 */

public class AtomSequencerCell extends AtomLinker {
    
    NeighborCell cell;       //cell currently occupied by this atom
    final AtomLinker nbrLink;  //linker used to arrange atom in sequence within the cells
    
    public AtomSequencerCell(Atom a) {
        super(a);
        nbrLink = new AtomLinker(a);
    }
    
    /**
     * Calls superclass method to remove atom from child list of parent,
     * and then removes atom from cell's occupant list.
     */
    public void remove() {
        super.remove();
        if (cell != null) {
            cell.occupants().remove(nbrLink);
            cell = null;
        }
    }
    
    /**
     * @return Returns the cell.
     */
    public NeighborCell getCell() {
        return cell;
    }
    /**
     * Singleton factory suitable to passing to the Atom constructor to specify
     * that atom sequencers should be this class.
     */
    public static final AtomSequencerFactory FACTORY = new AtomSequencerFactory() {
        public AtomLinker makeSequencer(Atom atom) {return new AtomSequencerCell(atom);}
    };
}