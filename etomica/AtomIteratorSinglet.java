package etomica;

/**
 * Iterator that expires after returning a single atom, which is
 * specified by a call to the reset(Atom) method.
 * A second atom may be specified using the reset(Atom, Atom) method. 
 * Then if the initiation flag indicates SKIP_FIRST, the second atom
 * is the only one returned; if the flag indicates INCLUDE_FIRST, then
 * the first atom is the only one returned.
 *
 * @author David Kofke
 */
public class AtomIteratorSinglet extends AtomIteratorAbstract { //could equally well write to implement AtomIterator
    
    private Atom atom1, atom2, next;
    
    public AtomIteratorSinglet() {this(AtomIterator.INCLUDE_FIRST);}
    public AtomIteratorSinglet(AtomIterator.Initiation init) {
        super(init);
    }
    
    /**
     * Returns the atom given in the last call to reset(Atom).
     */
    public Atom defaultFirstAtom() {return atom1;}
    
    /**
     * Returns the atom given in the last call to reset(Atom).
     */
    public Atom defaultLastAtom() {
        return (initiation == AtomIterator.INCLUDE_FIRST) ? atom1 : atom2;}
    
    public boolean isOrdered(Atom a1, Atom a2) {
        if(a1 == null || a2 == null) return false;
        else if(a1 == a2 && (a1 == atom1 || a1 == atom2)) return true;
        else if(a1 == atom1 && a2 == atom2) return true;
        else return false;
    }
    
    /**
     * Returns true if the given atom is the atom passed to the last call to reset(Atom),
     * or is one of the atoms passed to the last call to reset(Atom, Atom).
     */
    public boolean contains(Atom a) {return (a == atom1 || a == atom2) && a != null;}
    
    /**
     * Resets iterator to return atom specified by previous call to reset(Atom), or
     * to the second atom of the last call to reset(Atom, Atom) if SKIP_FIRST.
     */
    public Atom reset() {
        next = (initiation == AtomIterator.INCLUDE_FIRST) ? atom1 : atom2;
        hasNext = (next != null); 
        return atom;
    }
    
    /**
     * Sets the given atom as the one returned by the iterator if INCLUDE_FIRST.
     */
    public Atom reset(Atom a) {
        atom1 = a;
        atom2 = null;
        return reset();
    }
    
    /**
     * Same as call to reset(first).  Second argument is ignored.
     */
    public Atom reset(Atom first, Atom last) {
        atom1 = first;
        atom2 = last;
        return reset();
    }
    
    public Atom next() {hasNext = false; return next;}
    
    public void allAtoms(AtomAction act) {
        reset();
        act.actionPerformed(next);
    }
}//end of AtomIterator.Singlet
        
