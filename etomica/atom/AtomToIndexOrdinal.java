package etomica.atom;

import java.io.Serializable;


/**
 * Defines the index as the Atom's ordinal.
 * @author andrew
 */
public class AtomToIndexOrdinal implements AtomToIndex, Serializable {

    /**
     * @throws NullPointerException if the atom is null.
     */
    public int getIndex(Atom atom) {
        return atom.node.getOrdinal()-1;
    }
    
    private static final long serialVersionUID = 1L;

}
