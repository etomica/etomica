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
        public Atom.Iterator makeAtomIteratorUp() {return new Atom.Iterator.Up();}
        public Atom.Iterator makeAtomIteratorDown() {return new Atom.Iterator.Down();}
        public AtomPair.Iterator.A makeAtomPairIteratorUp() {return new AtomPair.Iterator.Up(phase);}
        public AtomPair.Iterator.A makeAtomPairIteratorDown() {return new AtomPair.Iterator.Down(phase);}
        public AtomPair.Iterator makeAtomPairIteratorFull() {return new AtomPair.Iterator.Full(phase);}
//        public AtomPair.Iterator makeAtomPairIteratorHalfUp() {return new AtomPair.Iterator.HalfUp(phase);}
//        public AtomPair.Iterator makeAtomPairIteratorHalfDown() {return new AtomPair.Iterator.HalfDown(phase);}
        public void clear() {}
    //    public abstract AtomPair.Iterator makeUpAtomPairIterator();
    //    public abstract AtomPair.Iterator makeDownAtomPairIterator();
    }
/*
    public final simulate.AtomPair.Iterator.A makePairIteratorFull(Atom iF, Atom iL, Atom oF, Atom oL) {return parentSimulation.space.makePairIteratorFull(boundary,iF,iL,oF,oL);}
    public final simulate.AtomPair.Iterator.A makePairIteratorHalf(Atom iL, Atom oF, Atom oL) {return parentSimulation.space.makePairIteratorHalf(boundary,iL,oF,oL);}
    public final simulate.AtomPair.Iterator.A makePairIteratorFull() {return parentSimulation.space.makePairIteratorFull(boundary);}
    public final simulate.AtomPair.Iterator.A makePairIteratorHalf() {return parentSimulation.space.makePairIteratorHalf(boundary);}
    public final simulate.AtomPair.Iterator.A makeUpIterator() {return parentSimulation.space.makeUpNeighborIterator(boundary);}
    public final simulate.AtomPair.Iterator.A makeDownIterator() {return parentSimulation.space.makeDownNeighborIterator(boundary);}
*/