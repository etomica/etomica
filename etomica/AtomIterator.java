package etomica;
 
/**
* Base interface for atom iterators.
* Atom iterators yield a sequence of atoms in successive calls to the next() method.
* When all atoms have been returned, hasNext() returns false.
* Iterators are often defined to progress "Up" or "Down" the set of atoms.
* "Up" and "Down" are arbitrary designations, except that the iterators guarantee
* that if atom 1 is "up list" of atom 2, then atom 2 is "down list" of atom 1.
* "Up" and "Down" relation between any atoms may change during the course of the 
* simulation, but at any instant the order is consistent and reproducible.
* "Neighbor" iterators yield only atoms that are considered to be "neighbors" of
* a specified atom.  The definition of "neighbor" depends on the iterator.  "Up neighbors"
* are those neighbors uplist of the atom; likewise with "Down neighbors".
*
* @see IteratorFactory
*/
public abstract class AtomIterator implements java.io.Serializable {        

    protected boolean hasNext = false;
    
    /**
     * @return true if the iterator will return another atom with a subsequent call to next(), false otherwise
     */
    public boolean hasNext() {return hasNext;}

    /**
     * @return the next atom in the list
     */
    public abstract Atom next();

    public abstract void reset(IteratorDirective id);
    /**
        * Resets the iterator in reference to the given atom.
        * Exactly how the given atom affects the reset depends on the particular iterator.
        * 
        * @param a
        */
    public abstract void reset(Atom a);

    /**
        * Resets the iterator, so that it is ready to go through its list again.
        */
    public abstract void reset();

    /**
        * Performs the given Action on each atom in the list in sequence.
        * 
        * @param act
        * @see Atom.Action
        */
    public abstract void allAtoms(AtomAction act);
            

  //************** END OF METHODS FOR AtomIterator *****************//
  
  
  //************** BEGIN DECLARATIONS OF INNER CLASSES *************//
  
   /**
    * Generic iterator that permits addition and removal of atoms.
    */
    public static final class List extends AtomIterator {
        private Atom.Linker first, next;
        public List() {hasNext = false;}
        public void reset() {
            next = first;
            hasNext = (next != null);
        }
        /**
         * reset with atom argument not properly implemented yet
         */
        public void reset(Atom a) {reset();}
        public Atom next() { //does not check that next is non-null
            Atom atom = next.atom();
            next = next.next();
            hasNext = (next != null);
            return atom;
        }
        public void allAtoms(AtomAction act) {
            for(Atom.Linker link=first; link!=null; link=link.next()) {
                act.actionPerformed(link.atom());
            }
        }
        /**
        * Adds an atom to the set of atoms given by this iterator
        */
        public void addAtom(Atom a) {
            first = new Atom.Linker(a,first);
        }
        //will someday add a removeAtom method
           
    }//end of Atom.Iterator.List
     
   /**
    * Iterator that expires after returning a single atom
    */
    public static final class Singlet extends AtomIterator {
        private Atom atom;
        public Singlet() {hasNext = false;}
        public Singlet(Atom a) {reset(a);}
        public void reset() {hasNext = (atom != null);}
        public void reset(Atom a) {atom = a; reset();}
        public Atom next() {hasNext = false; return atom;}
        public void allAtoms(AtomAction act) {act.actionPerformed(atom);}
    }
        
    /**
     * Iterator that progresses up a list of atoms.
     * Order of atoms is that given by the linked list of atoms, which changes
     * only if atoms are added or removed to/from the phase.
     */
    public static class Up extends AtomIterator  {
        protected Atom atom, nextAtom;
        private Phase phase;
        /**
         * Sets phase but does not reset iterator.  
         * Initial state is hasNext = false.
         */
        public Up(Phase p) {phase = p; hasNext = false;}
        /**
         * Sets phase and resets iterator using given atom.  
         * Initial state is hasNext = true if atom is not null.
         */
        public Up(Phase p, Atom a) {phase = p; reset(a);}
        public Phase phase() {return phase;}
        /**
         * Sets the iterator so the next atom is the one given (which may be null)
         */
        public void reset(Atom a) {
            atom = a;
            hasNext = (a != null);
        }
        /**
         * Sets the iterator so that the next atom is the first atom of the phase
         */
        public void reset() {reset(phase.firstAtom());}
        public Atom next() {
            nextAtom = atom;
            atom = atom.nextAtom();
            hasNext = (atom != null);
            return nextAtom;
        }
        /**
         * Performs the given action on all atoms in the phase.  Unaffected by any prior calls to reset.
         */
        public void allAtoms(AtomAction act) {
            for(Atom a=phase.firstAtom(); a!=null; a=a.nextAtom()) {act.actionPerformed(a);}
        }
    } //end of Atom.Iterator.Up
        
    /**
        * Iterates over all neighbors uplist from given atom.
        * This iterator defines <i>all</i> atoms in phase as neighbors, so it
        * simply extends the Up iterator, modifying the reset method to start with
        * the atom following the given atom, rather than the given atom itself.
        * Also, the no-argument reset performs a reset using the current atom, rather than 
        * resetting to neighbor of first atom in phase.
        */
    public static final class UpNeighbor extends Up {
        private Atom first;
        public UpNeighbor(Phase p) {super(p);}
        public UpNeighbor(Phase p, Atom a) {super(p,a);} 
        /**
            * Resets iterator so that the next atom is the one just upList of the given atom.
            */
        public void reset(Atom a) {
            atom = a;
            if(a == null) {hasNext = false; return;}
            first = a.nextAtom();
            super.reset(first);
        }
        /**
            * Resets iterator to the condition it was in after the last call to reset(Atom a).
            * This will be hasNext = false if reset(Atom a) was not called previously.
            */
        public void reset() {super.reset(first);}
        /**
            * Performs the given action on all atoms uplist of the one indicated in the last call to reset(Atom).  
            * If reset has not been called before, performs no action.
            */
        public void allAtoms(AtomAction act) {
            if(first == null) return;
            for(Atom a=first; a!=null; a=a.nextAtom()) {act.actionPerformed(a);}
        }
    }
        
    /**
        * Iterator that progresses down a list of atoms.
        * Order of atoms is that given by the linked list of atoms, which changes
        * only if atoms are added or removed from the phase.
        */
    public static class Down extends AtomIterator  {
        protected Atom atom;
        private Phase phase;
        /**
            * Sets phase but does not reset iterator.  
            * Initial state is hasNext = false.
            */
        public Down(Phase p) {phase = p; hasNext = false;}
        /**
            * Sets phase and resets iterator using given atom.  
            * Initial state is hasNext = true if atom is not null.
            */
        public Down(Phase p, Atom a) {phase = p; reset(a);}
        public Phase phase() {return phase;}
        public void reset(Atom a) {
            atom = a;
            hasNext = (a != null);
        }
        /**
            * Resets iterator to the first atom of the list.
            * Iterator will return only this atom and then expire, since there is nothing downlist of it
            */
        public void reset() {reset(phase.firstAtom());}
        public Atom next() {
            Atom nextAtom = atom;
            atom = atom.previousAtom();
            if(atom == null) {hasNext = false;}
            return nextAtom;
        }
        /**
            * Performs the given action on all atoms in the phase, starting from the last to the first.
            */
        public void allAtoms(AtomAction act) {
            for(Atom a=phase.lastAtom(); a!=null; a=a.previousAtom()) {act.actionPerformed(a);}
        }
    } //end of Atom.Iterator.Down
        
    /**
        * Iterates over all neighbors downlist from given atom.
        * This iterator defines <i>all</i> atoms in phase as neighbors, so it
        * simply extends the Down iterator, modifying the reset method to start with
        * the atom preceding the given atom, rather than the given atom itself.
        * Also, the no-argument reset performs a reset using the current atom, rather than 
        * resetting to neighbor of first atom in phase.
        */
    public static final class DownNeighbor extends Down {
        private Atom first;
        public DownNeighbor(Phase p) {super(p);}
        public DownNeighbor(Phase p, Atom a) {super(p,a);}
        /**
            * Resets iterator so that the next atom is the one just downList of the given atom
            */
        public void reset(Atom a) {
            atom = a;
            if(a == null) {hasNext = false; return;}
            first = a.previousAtom();
            super.reset(first);
        }
        /**
            * Resets iterator to the condition it was in after the last call to reset(Atom a)
            * This will be hasNext = false if reset(Atom a) was not called previously.
            */
        public void reset() {super.reset(first);}
        /**
            * Performs the given action on all atoms uplist of the one indicated in the last call to reset(Atom).  
            * If reset has not been called before, performs no action.
            */
        public void allAtoms(AtomAction act) {
            if(first == null) return;
            for(Atom a=first; a!=null; a=a.previousAtom()) {act.actionPerformed(a);}
        }
    }//end of DownNeighbor
}//end of AtomIterator
    
