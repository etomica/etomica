package etomica.atom;

import java.io.Serializable;

import etomica.atom.iterator.AtomIteratorArraySequence.AtomToIndex;

/**
 * Defines the index as the Atom's ordinal.
 * @author andrew
 */
public class AtomToIndexOrdinal implements AtomToIndex, Serializable {

    public int getIndex(Atom atom) {
        return (atom != null) ? atom.node.getOrdinal() : -1;
    }
    
    private static final long serialVersionUID = 1L;

}
