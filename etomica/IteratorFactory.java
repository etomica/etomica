package etomica;

/**
 * Class for construction of Atom iterators.
 * Default iterators use the atom linked lists in a phase.
 * This class can be extended to define iterators that use, for example, neighbor lists.
 * "Neighbor" iterators yield only atoms that are considered to be "neighbors" of
 * a specified atom.  The definition of "neighbor" depends on the iterator.  "Up neighbors"
 * are those neighbors uplist of the atom; likewise with "Down neighbors".
 *
 * @author David Kofke
 */

public interface IteratorFactory {

    /**
     * creates a pair iterator which loops over all pairs in a neighbor list
     * between two groups
     * @return the pair iterator
     */
    public AtomsetIterator makeIntergroupNbrPairIterator();
    
    /**
     * creates a pair iterator which loops over all pairs in a neighbor list
     * within one group
     * @return the pair iterator
     */
    public AtomsetIterator makeIntragroupNbrPairIterator();
    
    /**
     * Returns an atom sequencer that does not base ordering on any type
     * of neighbor scheme.  Sequencer is assigned by atom factory when
     * creating a new atom.
     */
    public AtomSequencer makeSimpleSequencer(Atom atom);
    
    /**
     * Returns an atom sequencer that bases ordering on a neighbor scheme.
     * The nature of the neighbor scheme depends on the way neighboring is
     * set up by the IteratorFactory products. Sequencer is assigned by 
     * atom factory when creating a new atom.
     */
    public AtomSequencer makeNeighborSequencer(Atom atom);
    
    public Class simpleSequencerClass();
    
    public Class neighborSequencerClass();
    
    public AtomSequencer.Factory simpleSequencerFactory();
    
    public AtomSequencer.Factory neighborSequencerFactory();
   
}
