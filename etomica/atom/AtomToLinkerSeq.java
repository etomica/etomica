package etomica.atom;

import java.io.Serializable;

import etomica.atom.iterator.AtomIteratorSequence.AtomToLinker;

/**
 * Defines the linker as the Atom's sequencer.
 * @author andrew
 */
public class AtomToLinkerSeq implements AtomToLinker, Serializable {

    public AtomLinker getLinker(Atom atom) {
        return (atom != null) ? atom.seq : null;
    }

}
