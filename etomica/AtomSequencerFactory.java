/*
 * History
 * Created on Dec 2, 2004 by kofke
 */
package etomica;

/**
 * AtomSequencer is used to structure all the atoms in a phase into
 * a well defined order.  A single instance of this class is held by
 * the seq field of each atom, and it is the primary point of reference
 * for structuring lists of child atoms in each atom group.  Most of 
 * the iterators that loop through the atoms in the phase use the list
 * order set up using the sequencer.
 *
 * @author David Kofke
 * @version 02.03.09
 */
public interface AtomSequencerFactory {
    
    public AtomLinker makeSequencer(Atom atom);
    
    public static final AtomSequencerFactory SIMPLE = new AtomSequencerFactory() {
        public AtomLinker makeSequencer(Atom atom) {
            return new AtomLinker(atom);
        }
    };

}