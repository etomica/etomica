package etomica;
 
/**
* Parent class for atom iterators.
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
* @author David Kofke
*/
public abstract class AtomIterator implements java.io.Serializable {        

    protected boolean hasNext;
    
    /**
     * Iterator is constructed not ready for iteration.  Must call a reset
     * method before use.  hasNext returns false until then.
     */
    public AtomIterator() {
        hasNext = false;
    }
    
    /**
     * @return true if the iterator will return another atom with a subsequent 
     * call to next(), false otherwise.
     */
    public boolean hasNext() {return hasNext;}

    public void reset(IteratorDirective id) {
        IteratorDirective.Bounds bounds = id.bounds();
        if(bounds == IteratorDirective.ALL)             reset();
        else if(bounds == IteratorDirective.FIRST)      reset(id.firstAtom());
        else if(bounds == IteratorDirective.FIRST_LAST) reset(id.firstAtom(), id.lastAtom());
        else hasNext = false;
    }

    /**
     * @return the next atom in the list
     */
    public abstract Atom next();
    
    /**
     * Resets the iterator, so that it is ready to go through its list again, beginning
     * with its natural first iterate (the identity of which depends on the iterator).
     *
     * @return the atom that will be returned with the first call to next()
     */
    public abstract Atom reset();

    /**
     * Resets the iterator in reference to the given atom.  If the atom is among the iterates
     * of this iterator, then it will be the next one returned on a call to next().  If
     * it is not one of the iterates, this is equivalent to a call to reset().
     * 
     * @param first  the nominal first atom returned by the iterator
     * @return       the atom that will be returned with the next subsequent call to next()
     */
    public abstract Atom reset(Atom first);

    /**
     * Resets iterator in reference to the given atoms.  Initiation of the iteration is
     * as given by a call to reset(first), and termination of th iteration is done in reference
     * to the second given atom.  If this atom is among the iterates of this iterator, iteration
     * will terminate when this atom is reached and returned by next(); if it is not one of
     * the iterates, iteration terminates when the nature end of the iteration is reached.
     *
     * @param first the nominal first atom returned by the iterator
     * @param last  the nominal last atom returned by the iterator
     * @return      the atom that will be returned with the first call to next()
     */
    public abstract Atom reset(Atom first, Atom last);
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
        private AtomLinker first, last, next;
        private Atom terminator;
        public List() {super();}
        
        /**
         * Sets to iterate over all atoms in list.
         */
        public Atom reset() {
            next = first;
            hasNext = (next != null);
            if(!hasNext) return null;
            terminator = last.atom();
            return next.atom();
        }
        /**
         * Sets to begin iterating from given atom to end of list.  If atom is not in list
         * set to iterate from beginning of list.
         */
        public Atom reset(Atom firstAtom) {
            //look for firstAtom in list
            next = first;
            while(next != null && next.atom() != firstAtom) next = next.next();
            if(next == null) next = first; //firstAtom is not in list
            
            hasNext = (next != null);
            if(!hasNext) return null;
            terminator = last.atom();
            return next.atom();
        }
        /**
         * Sets to begin iterating between given atoms.
         */
        public Atom reset(Atom firstAtom, Atom lastAtom) {
            Atom nextAtom = reset(firstAtom);
            terminator = lastAtom;
            return nextAtom;
       /*     if(first == null) {hasNext = false; return null;} //empty list
            next = first;
            while(next.atom() != firstAtom) {
                if(next == last || next.atom() == lastAtom) {//reached end of list, or encountered lastAtom before firstAtom
                    hasNext = false; return null;
                }
                next = next.next();
            }
            AtomLinker next2 = next;
            while(next2 != null && next2.atom() != lastAtom) next2 = next2.next();
            if(next2 == null) {hasNext = false; return null;} //lastAtom is not in list
            terminator = lastAtom;
            hasNext = true;
            return next.atom(); */
        }
        
        public Atom next() { //does not check that next is non-null
            Atom atom = next.atom();
            next = next.next();
            hasNext = (atom != terminator);
            return atom;
        }
        /**
         * Performs action on all atoms in list.  Not affected by iteratorDirective given in reset.
         */
        public void allAtoms(AtomAction act) {
            for(AtomLinker link=first; link!=null; link=link.next()) {
                act.actionPerformed(link.atom());
            }
        }
        /**
        * Adds an atom to the set of atoms given by this iterator
        */
        public void addAtom(Atom a) {
            if(a == null) return;
            AtomLinker newLink = new AtomLinker(a);
            if(last != null) last.setNext(newLink);
            last = newLink;
            if(first == null) first = newLink;
        }
        //will someday add a removeAtom method
           
    }//end of AtomIterator.List
     
   /**
    * Iterator that expires after returning a single atom
    */
    public static final class Singlet extends AtomIterator {
        private Atom atom;
        public Singlet() {super();}
        public Singlet(Atom a) {reset(a);}
        /**
         * Resets iterator to return atom specified by previous call to reset(Atom).
         */
        public Atom reset() {hasNext = (atom != null); return atom;}
        /**
         * Sets the given atom as the one returned by the iterator.
         */
        public Atom reset(Atom a) {atom = a; return reset();}
        /**
         * Same as call to reset(first).  Second argument is ignored.
         */
        public Atom reset(Atom first, Atom last) {return reset(first);}
        public Atom next() {hasNext = false; return atom;}
        public void allAtoms(AtomAction act) {act.actionPerformed(atom);}
    }//end of AtomIterator.Singlet
        
    /**
     * Iterator that progresses up a list of atoms.
     * Order of atoms is that given by the linked list of atoms, which changes
     * only if atoms are added or removed to/from the phase.
     */
    public static class Up extends AtomIterator  {
        protected Atom atom, nextAtom, terminator;
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
    } //end of AtomIterator.Up
        
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
        
    /**
     * Iterator that progresses down a list of atoms.
     * Order of atoms is that given by the linked list of atoms, which changes
     * only if atoms are added or removed from the phase.
     */
    public static class Down extends AtomIterator  {
        protected Atom atom, terminator;
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
        public Atom reset(Atom a) {
            atom = a;
            terminator = null;
            hasNext = (a != null);
            return atom;
        }
        /**
         * Sets the iterator to loop through the given atoms, inclusive.
         * Does not check that last atom will be encountered.  If not, iteration
         * will terminate when last iterable atom is reached.
         */
        public Atom reset(Atom first, Atom last) {
            reset(first);
            terminator = last;
            return atom;
        }
        /**
         * Resets iterator to the first atom of the list.
         * Iterator will return only this atom and then expire, since there is nothing downlist of it
         */
        public Atom reset() {return reset(phase.firstAtom());}
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
    } //end of AtomIterator.Down
        
    /**
        * Iterates over all neighbors downlist from given atom.
        * This iterator defines <i>all</i> atoms in phase as neighbors, so it
        * simply extends the Down iterator, modifying the reset method to start with
        * the atom preceding the given atom, rather than the given atom itself.
        * Also, the no-argument reset performs a reset using the current atom, rather than 
        * resetting to neighbor of first atom in phase.
        */
 /*   public static final class DownNeighbor extends Down {
        private Atom first;
        public DownNeighbor(Phase p) {super(p);}
        public DownNeighbor(Phase p, Atom a) {super(p,a);}
        /**
            * Resets iterator so that the next atom is the one just downList of the given atom
            * /
        public void reset(Atom a) {
            atom = a;
            if(a == null) {hasNext = false; return;}
            first = a.previousAtom();
            super.reset(first);
        }
        /**
         * Resets iterator to the condition it was in after the last call to reset(Atom a)
         * This will be hasNext = false if reset(Atom a) was not called previously.
         * /
        public void reset() {super.reset(first);}
        /**
         * Performs the given action on all atoms uplist of the one indicated in the last call to reset(Atom).  
         * If reset has not been called before, performs no action.
         * /
        public void allAtoms(AtomAction act) {
            if(first == null) return;
            for(Atom a=first; a!=null; a=a.previousAtom()) {act.actionPerformed(a);}
        }
    }//end of AtomIterator.DownNeighbor
    */
}//end of AtomIterator
    
