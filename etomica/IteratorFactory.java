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
     * Returns iterator that loops over all atoms in a group without
     * reference to another atom and using any convenient order.  Does not
     * necessarily adhere to up/down ordering.  This iterator is required
     * if looping over atoms while performing operations on them that might
     * cause the sequencing to be shuffled during iteration.  In general
     * it is also more simply implemented and presumably faster than the
     * analogous sequential iterator (given by makeGroupIteratorSequential).
     */
    //public AtomIterator makeGroupIteratorSimple();
    
    /**
     * Returns iterator that loops over all atoms in a group, without
     * reference to another atom.  Adheres to atom sequence, so that
     * instructions for up or down iteration are followed using ordering
     * consistent with that followed by neighbor iterators (e.g., IntergroupIterator).
     */
    public AtomIterator makeGroupIteratorSequential();
    
    /**
     * Returns iterator that loops over all atoms in a group, with reference
     * to an atom in another group.  The atoms returned by the iterator
     * would then be "neighbors" of the given atom.  Definition of "neighbor"
     * depends on iteration scheme defined by the IteratorFactory.  Atoms
     * are returned in an order that is consistent with that given by
     * the iterator from makeGroupIteratorSequential (except that this
     * neighbor iterator may return fewer atoms).
     */
    public AtomIterator makeIntergroupNbrIterator();

    /**
     * Returns iterator that loops over all atoms in a group, with reference
     * to one of the atoms in the group.  The atoms returned by the iterator
     * would then be "neighbors" of the given atom.  Definition of "neighbor"
     * depends on iteration scheme defined by the IteratorFactory.
     */
    public AtomIterator makeIntragroupNbrIterator();

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
