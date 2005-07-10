package etomica.atom;

import etomica.Atom;
import etomica.Debug;
import etomica.Phase;

/**
 * AtomSource that returns a completely random leaf atom.
 */
public class AtomSourceRandomLeaf implements AtomSource, java.io.Serializable {

    public void setPhase(Phase p) {
        list = p.getSpeciesMaster().atomList;
    }
    
    /**
     * returns a random atom from the phase's leaf atom list
     */
    public Atom getAtom() {
        if (Debug.ON && list== null) throw new IllegalStateException("must set the phase before calling getAtom");
        return list.getRandom();
    }
    
    protected AtomList list = null;
}
