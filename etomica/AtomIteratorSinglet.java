package etomica;

/**
 * Iterator that expires after returning a single atom.
 *
 * @author David Kofke
 */
public final class AtomIteratorSinglet extends AtomIterator {
    
    private Atom atom;
    
    public AtomIteratorSinglet() {super();}
    public AtomIteratorSinglet(Atom a) {reset(a);}
    
    /**
     * Returns true if the given atom is the atom passed to that last call to reset(Atom).
     */
    public boolean contains(Atom a) {return a == atom;}
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
        
