package etomica;

/**
 * Iterator that progresses up a list of atoms.
 * Order of atoms is that given by the linked list of atoms, which changes
 * only if atoms are added or removed to/from the phase.
 */
public class AtomIteratorUp extends AtomIterator  {
    protected Atom atom, nextAtom, terminator;
    private Phase phase;
    /**
     * Sets phase but does not reset iterator.  
     * Initial state is hasNext = false.
     */
    public AtomIteratorUp(Phase p) {phase = p; hasNext = false;}
    /**
     * Sets phase and resets iterator using given atom.  
     * Initial state is hasNext = true if atom is not null.
     */
    public AtomIteratorUp(Phase p, Atom a) {phase = p; reset(a);}
    public Phase phase() {return phase;}
        
    //revisit this
    public boolean contains(Atom atom) {return true;}

    /**
     * Sets the iterator so the next atom is the one given (which may be null).
     * Iteration proceeds through the last iterable atom.
     */
    public Atom reset(Atom a) {
        atom = a;
        terminator = null;
        hasNext = (a != null);
        return atom;
    }
    /**
     * Sets the iterator so that the next atom is the first atom of the phase
     */
    public Atom reset() {return reset(phase.firstAtom());}
    /**
     * Sets the iterator to loop through the given atoms, inclusive.
     * Does not check that last atom will be encountered.  If it isn't, iteration
     * will terminate when last iterable atom is reached.
     */
    public Atom reset(Atom first, Atom last) {
        reset(first);
        terminator = last;
        return atom;
    }
    public Atom next() {
        nextAtom = atom;
        atom = atom.nextAtom();
        hasNext = (atom != null && nextAtom != terminator);
        return nextAtom;
    }
    /**
     * Performs the given action on all atoms in the phase.  Unaffected by any prior calls to reset.
     */
    public void allAtoms(AtomAction act) {
        for(Atom a=phase.firstAtom(); a!=null; a=a.nextAtom()) {act.actionPerformed(a);}
    }
} //end of AtomIteratorUp

        
    /**
        * Iterates over all neighbors uplist from given atom.
        * This iterator defines <i>all</i> atoms in phase as neighbors, so it
        * simply extends the Up iterator, modifying the reset method to start with
        * the atom following the given atom, rather than the given atom itself.
        * Also, the no-argument reset performs a reset using the current atom, rather than 
        * resetting to neighbor of first atom in phase.
        */
/*    public static final class UpNeighbor extends Up {
        private Atom first;
        public UpNeighbor(Phase p) {super(p);}
        public UpNeighbor(Phase p, Atom a) {super(p,a);} 
        /**
         * Resets iterator so that the next atom is the one just upList of the given atom.
         * /
        public void reset(Atom a) {
            atom = a;
            if(a == null) {hasNext = false; return;}
            first = a.nextAtom();
            super.reset(first);
        }
        /**
         * Resets iterator to the condition it was in after the last call to reset(Atom a).
         * This will be hasNext = false if reset(Atom a) was not called previously.
         * /
        public void reset() {super.reset(first);}
        /**
         * Performs the given action on all atoms uplist of the one indicated in the last call to reset(Atom).  
         * If reset has not been called before, performs no action.
         * /
        public void allAtoms(AtomAction act) {
            if(first == null) return;
            for(Atom a=first; a!=null; a=a.nextAtom()) {act.actionPerformed(a);}
        }
    }//end of AtomIterator.UpNeighbor
*/