package etomica;

/**
 * Generic iterator that permits addition and removal of atoms.
 *
 * @author David Kofke
 */
public final class AtomIteratorList extends AtomIterator {
    
    private AtomLinker first, last, next;
    private Atom terminator;
    
    public AtomIteratorList() {super();}
        
    public boolean contains(Atom atom) {
        AtomLinker link = first;
        while(link != null && link.atom() != atom) link = link.next();
        return (link != null && link.atom() == atom);
    }

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
        hasNext = (atom != terminator && next != null);
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
