package etomica;        
/**
 * Iterator that progresses down a list of atoms.
 * Order of atoms is that given by the linked list of atoms, which changes
 * only if atoms are added or removed from the phase.
 *
 * @author David Kofke
 */
public class AtomIteratorDown extends AtomIterator  {
    
    protected Atom atom, terminator;
    private Phase phase;
    /**
     * Sets phase but does not reset iterator.  
     * Initial state is hasNext = false.
     */
    public AtomIteratorDown(Phase p) {
        phase = p; 
        hasNext = false;
    }
    /**
     * Sets phase and resets iterator using given atom.  
     * Initial state is hasNext = true if atom is not null.
     */
    public AtomIteratorDown(Phase p, Atom a) {
        phase = p; 
        reset(a);
    }
    
    public Phase phase() {return phase;}
        
    //revisit this
    public boolean contains(Atom atom) {return true;}
        
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
} //end of AtomIteratorDown
        
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
