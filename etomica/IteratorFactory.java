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
        public AtomIterator makeAtomIterator() {return phase.new Iterator();}
/*        public AtomIterator makeAtomIteratorUp() {return new AtomIteratorUp(phase);}
        public AtomIterator makeAtomIteratorDown() {return new AtomIteratorDown(phase);}
        public AtomIterator makeAtomIteratorUpDown() {return new AtomIteratorUpDown(phase);}
        public AtomIterator makeAtomIterator(AtomIterator.Initiation init) {return new AtomIteratorUp(phase, init);}
        public AtomIterator makeAtomIteratorUp(AtomIterator.Initiation init) {return new AtomIteratorUp(phase, init);}
        public AtomIterator makeAtomIteratorDown(AtomIterator.Initiation init) {return new AtomIteratorDown(phase, init);}
        public AtomIterator makeAtomIteratorUpDown(AtomIterator.Initiation init) {return new AtomIteratorUpDown(phase, init);}
        public AtomIterator makeAtomIteratorUpNeighbor() {return new AtomIteratorUp(phase, AtomIterator.SKIP_FIRST);}
        public AtomIterator makeAtomIteratorDownNeighbor() {return new AtomIteratorDown(phase, AtomIterator.SKIP_FIRST);}
  */
        public void addMolecule(Molecule m) {}
        public void deleteMolecule(Molecule m) {}
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
