package etomica;

/**
 * Class for construction of Atom iterators.
 * Default iterators use the atom linked lists in a phase.
 * This class can be extended to define iterators that use, for example, neighbor lists.
 * "Neighbor" iterators yield only atoms that are considered to be "neighbors" of
 * a specified atom.  The definition of "neighbor" depends on the iterator.  "Up neighbors"
 * are those neighbors uplist of the atom; likewise with "Down neighbors".
 * Each instance of a Phase has one instance of an Iterator for iterating over its atoms
 *
 * @author David Kofke
 * @see Atom.Iterator
 * @see AtomPair.Iterator
 */

    public class IteratorFactory implements java.io.Serializable {
        private Phase phase;
        public IteratorFactory(Phase p) {phase = p;}
        public AtomIterator makeAtomIterator() {return phase.makeAtomIterator();}

        public void moveNotify(Atom a) {}  //notification by integrator when an atom is moved
        public void reset() {}
        public Phase phase() {return phase;}
   
      /**
       * Interface for a class that can make an iterator factory
       */
        public interface Maker extends java.io.Serializable {
            public IteratorFactory makeIteratorFactory(Phase p);
        }
    }
