package etomica;

/**
 * Class for construction of Atom iterators.  Makes single-atom and pair-atom iterators.
 * Default iterators are those defined in Atom and AtomPair, which use the atom linked lists
 * This class can be extended to define iterators that use, for example, neighbor lists.
 * Each instance of a Phase has one instance of an Iterator for iterating over its atoms
 *
 * @see Atom.Iterator
 * @see AtomPair.Iterator
 */

    public class IteratorFactory implements java.io.Serializable {
        private Phase phase;
        public IteratorFactory(Phase p) {phase = p;}
        public AtomIterator makeAtomIterator() {return new AtomIteratorUp(phase);}
        public AtomIterator makeAtomIteratorUp() {return new AtomIteratorUp(phase);}
        public AtomIterator makeAtomIteratorDown() {return new AtomIteratorDown(phase);}
        public AtomIterator makeAtomIteratorUpNeighbor() {return new AtomIteratorUp(phase);}
        public AtomIterator makeAtomIteratorDownNeighbor() {return new AtomIteratorDown(phase);}
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
