package etomica;

/**
 * Iterator that progresses in both directions over a list of atoms.  Differs
 * from Up and Down iterators in that iterations
 * from an atom given via reset(Atom) proceed first up from it, then down. 
 * Also, iterations from one atom to another given via reset(Atom, Atom) proceeds
 * regardless of their order, as long as both are among the iterates.  Iteration
 * following reset() proceeds as in the Up iterator.
 * Order of atoms is that given by the linked list of atoms, which changes
 * only if atoms are added or removed to/from the phase.
 *
 * @author David Kofke
 */
public class AtomIteratorUpDown extends AtomIteratorAbstract  {

    private final Phase phase;
    private boolean upList, goDown;
    private Atom setAtom;
    
    /**
     * Sets phase but does not reset iterator.  
     * Initial state is hasNext = false.
     * Default initiation flag is INCLUDE_FIRST.
     */
    public AtomIteratorUpDown(Phase p) {
        this(p, AtomIterator.INCLUDE_FIRST);
    }
    public AtomIteratorUpDown(Phase p, AtomIterator.Initiation init) {
        super(init);
        phase = p; 
    }

    public Phase phase() {return phase;}
    
    /**
     * Returns first atom in phase.
     */
    public Atom defaultFirstAtom() {return phase.firstAtom();}
    
    /**
     * Returns last atom in phase.
     */
    public Atom defaultLastAtom() {return phase.lastAtom();}
    
    public boolean isOrdered(Atom atom1, Atom atom2) {
        return atom1.preceeds(atom2);
    }
    
    /**
     * Set to go up list from default first.
     */
    public Atom reset() {
        reset(defaultFirstAtom());
        goDown = false;
        return atom;
    }
    
    /**
     * Set to go up, then down from given atom.
     */
    public Atom reset(Atom a) {
        if(!contains(a)) {  //also ensures that a is not null
            hasNext = false;
            return null;
        }
        setAtom = a;
        atom = a;
        terminator = defaultLastAtom();
        if(initiation == SKIP_FIRST && atom != null) atom = atom.nextAtom();
        hasNext = (atom != null);
        upList = true;
        goDown = true;
        return atom;
    }
    
    private void resetDown() {
        atom = setAtom;
        terminator = defaultFirstAtom();
        if(atom != null) atom = atom.previousAtom(); //uplist would have handled first atom, so skip it
        hasNext = (atom != null);
        upList = false;
    }
    
    public Atom reset(Atom first, Atom last) {
        if(isOrdered(first,last)) { //also ensures non-null
            atom = first;
            terminator = last;
            if(initiation == SKIP_FIRST) atom = atom.nextAtom();
            hasNext = (atom != null);
            upList = true;
            goDown = false;
            return atom;
        }
        else if(isOrdered(last, first)) {
            atom = first;
            terminator = last;
            if(initiation == SKIP_FIRST) atom = atom.previousAtom();
            hasNext = (atom != null);
            upList = false;
            return atom;
        }
        else {
            hasNext = false;
            return null;
        }
    }

    public Atom next() {
        Atom nextAtom = atom;
        atom = upList ? atom.nextAtom() : atom.previousAtom();
        hasNext = nextAtom != terminator;
        if(!hasNext && upList && goDown) resetDown();
        return nextAtom;
    }
    /**
     * Performs the given action on all atoms in the phase.  Unaffected by any prior calls to reset.
     */
    public void allAtoms(final AtomAction act) {
        for(Atom a=phase.firstAtom(); a!=null; a=a.nextAtom()) {act.actionPerformed(a);}
    }
} //end of AtomIteratorUpDown


/*package etomica;


public final class AtomIteratorUpDown implements AtomIterator {
    
    private AtomIterator up, down;
rUpDown as in the Up itera, goDown;
    private boolean hasNext;
    
    public class AtomIteratorUpDown(AtomIterator up, AtomIterator down) {
        setIterators(up, down);
    }
    
    public void setIterators((AtomIterator up, AtomIterator down) {
        if(up != null) this.up = up;  //don't overwrite with a null
        if(down != null) this.down = down;
    }
    
    public boolean hasNext() {return hasNext;}
    
    public boolean contains(Atom a) {return up.contains(a) || down.contains(a);}
    
    public void reset(IteratorDirective id) {
        switch(id.atomCount()) {
            case 0:  reset(); 
                     break;
            case 1:  reset(id.atom1()); 
                     break;
            case 2:  reset(id.atom1(), id.atom2()); 
                     break;
            default: hasNext = false; 
                     break;
        }
    }
    
    public void reset() {
        up.reset();
        goDown = false;
    }
    
    */