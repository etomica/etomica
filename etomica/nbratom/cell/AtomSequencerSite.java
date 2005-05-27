/*
 * History
 * Created on Nov 23, 2004 by kofke
 */
package etomica.nbratom.cell;

import etomica.Atom;
import etomica.atom.AtomLinker;
import etomica.atom.AtomSequencerFactory;


/**
 * Sequencer used for atoms existing on a lattice.
 */

public class AtomSequencerSite extends AtomLinker {
    
    AtomSite site;       //site currently occupied by this atom
    
    public AtomSequencerSite(Atom a) {
        super(a);
    }
    
    /**
     * Calls superclass method to remove atom from child list of parent,
     * and then removes atom from cell's occupant list.
     */
    public void remove() {
        super.remove();
        if (site != null) {
            site.setAtom(null);
            site= null;
        }
    }
    
    /**
     * @return Returns the cell.
     */
    public AtomSite getSite() {
        return site;
    }
    /**
     * Singleton factory suitable to passing to the Atom constructor to specify
     * that atom sequencers should be this class.
     */
    public static final AtomSequencerFactory FACTORY = new AtomSequencerFactory() {
        public AtomLinker makeSequencer(Atom atom) {return new AtomSequencerSite(atom);}
    };
}