package simulate;

/**
 * Interface for construction of Atom iterators.  Makes single-atom and pair-atom iterators.
 * Default iterators are those defined in Atom and AtomPair, which use the atom linked lists
 * This class can be extended to define iterators that use, for example, neighbor lists
 * Each instance of a Phase has one instance of an Iterator for iterating over its atoms
 */

/*interface Iterator {
    public Atom.Iterator makeAtomIteratorUp();
    public Atom.Iterator makeAtomIteratorDown();
    public AtomPair.Iterator makeAtomPairIteratorFull(Phase p);
    public AtomPair.Iterator makeAtomPairIteratorHalf(Phase p);
*/
    public class Iterator {
        public static final boolean INTRA = true;
        public static final boolean INTER = false;
        public Phase phase;
        public Iterator(Phase p) {phase = p;}
        public Atom.Iterator makeAtomIteratorUp() {return new Atom.Iterator.Up(phase);}
        public Atom.Iterator makeAtomIteratorDown() {return new Atom.Iterator.Down(phase);}
        public AtomPair.Iterator.A makeAtomPairIteratorUp() {return new AtomPair.Iterator.Up(phase);}
        public AtomPair.Iterator.A makeAtomPairIteratorDown() {return new AtomPair.Iterator.Down(phase);}
        public AtomPair.Iterator makeAtomPairIteratorAll() {return new AtomPair.Iterator.All(phase);}
        public void addMolecule(Molecule m) {}
        public void deleteMolecule(Molecule m) {}
        public void moveNotify(Atom a) {}  //notification by integrator when an atom is moved
        public void reset() {}
   
      /**
       * Interface for a class that can make an iterator
       */
        interface Maker {
            public Iterator makeIterator(Phase p);
        }
    }
