package simulate;

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
        public static final boolean INTRA = true;
        public static final boolean INTER = false;
        private Phase phase;
        public IteratorFactory(Phase p) {phase = p;}
        public Atom.Iterator makeAtomIterator() {return makeAtomIteratorUp();}
        public Atom.Iterator makeAtomIteratorUp() {return new Atom.Iterator.Up(phase);}
        public Atom.Iterator makeAtomIteratorDown() {return new Atom.Iterator.Down(phase);}
        public Atom.Iterator makeAtomIteratorUpNeighbor() {return new Atom.Iterator.UpNeighbor(phase);}
        public Atom.Iterator makeAtomIteratorDownNeighbor() {return new Atom.Iterator.DownNeighbor(phase);}
        public AtomPair.Iterator makeAtomPairIteratorUp() {return new AtomPair.Iterator.Up(phase);}
        public AtomPair.Iterator makeAtomPairIteratorDown() {return new AtomPair.Iterator.Down(phase);}
        public AtomPair.Iterator makeAtomPairIteratorAll() {return new AtomPair.Iterator.All(phase);}
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
