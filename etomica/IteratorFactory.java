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
    public AtomsetIterator makeInterSpeciesPairIterator();
    
    /**
     * creates a pair iterator which loops over all pairs in a neighbor list
     * within one group
     * @return the pair iterator
     */
    public AtomsetIterator makeIntraSpeciesPairIterator();
    
    /**
     * Sequencer used for molecule-level atoms (those with a SpeciesAgent
     * as the parent). Special because if cell lists are used, they are
     * kept for these atoms.
      */
    public AtomSequencer.Factory moleculeSequencerFactory();
    
    /**
     * Sequencer used for atoms on which concrete potentials (non-group) act.
     * Special because if neighbor lists are used, they are kept for these atoms.
     * @return
     */
    public AtomSequencer.Factory interactionAtomSequencerFactory();

    /**
     * Sequencer used for molecule-level atoms on which concrete potentials 
     * (non-group) act.  Special because they must be able to handle both neighbor
     * and cell listing.
     * @return
     */
    public AtomSequencer.Factory interactionMoleculeSequencerFactory();

}
