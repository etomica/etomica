package etomica;

/**
 * Parent class for most atom iterators.
 *
 * @author David Kofke
 */
public abstract class AtomIteratorAbstract implements AtomIterator, java.io.Serializable {        

    protected boolean hasNext;
    protected Atom atom, nextAtom, terminator;
    protected boolean upListNow, doGoDown;
    protected boolean isNeighborIterator;
    protected Atom setAtom, terminator2;

    /**
     * Iterator is constructed not ready for iteration.  Must call a reset
     * method before use.  hasNext returns false until then.
     */
    public AtomIteratorAbstract() {
        hasNext = false;
    }
    
    /**
     * Natural first atom returned by this iterator.  First atom actually returned
     * may differ due to call to reset(Atom) or reset(Atom, Atom), or becuase 
     * initiation flag indicates skipping first atom.
     */
    public abstract Atom defaultFirstAtom();
    
    /**
     * Natural last atom returned by this iterator.  Last atom returned may differ
     * due to call to reset(Atom, Atom).
     */
    public abstract Atom defaultLastAtom();
    
    /**
     * Indicates if the given atom is among the iterates returned by this iterator.
     *
     * @return <code>true</code> if atom is one of the iterates.
     */
    public abstract boolean contains(Atom atom);

    /**
     * Returns true if atom1 preceeds atom2 in the up-direction sequence of iterates
     * returned by this iterator.  Also returns true if atom1 and atom2 are
     * the same.  Returns false if atom1 succeeds atom2 in sequence, or
     * if either atom1 or atom2 are not among the iterates in the list, or are null.
     */
    public abstract boolean isOrdered(Atom atom1, Atom atom2);
    
    /**
     * @return the next atom in the list
     */
    public abstract Atom next();
    
    protected abstract Atom firstUpNeighbor(Atom a);
    protected abstract Atom firstDownNeighbor(Atom a);
    
    /**
     * Performs the given Action on each atom in the list in sequence.
     * 
     * @param act
     * @see Atom.Action
     */
    public abstract void allAtoms(AtomAction act);
    
    /**
     * @return true if the iterator will return another atom with a subsequent 
     * call to next(), false otherwise.
     */
    public boolean hasNext() {return hasNext;}
    
    public void setAsNeighbor(boolean b) {isNeighborIterator = b;}

    /**
     * Returns the first neighbor of the given atom, considering
     * the value of the direction field.
     */
    public Atom firstNeighbor(Atom a) {
        Atom neighbor = null;
        if(upListNow) neighbor = firstUpNeighbor(a);
        if(neighbor == null && doGoDown) {
            upListNow = false;
            neighbor = firstDownNeighbor(a);
        }
        return neighbor;
    }
    
    public Atom reset(IteratorDirective id) {
        upListNow = id.direction().doUp();
        doGoDown = id.direction().doDown();
        switch(id.atomCount()) {
            case 0:  return reset(); 
            case 1:  return reset(id.atom1()); 
            case 2:  return reset(id.atom1(), id.atom2()); 
            default: hasNext = false; 
            return null;
        }
    }
    
    /**
     * Resets the iterator, so that it is ready to go through its list again, beginning
     * with its natural first iterate (the identity of which depends on the iterator).
     *
     * @return the atom that will be returned with the first call to next() 
     */
    public Atom reset() {
        if(isNeighborIterator) atom = null;
        else if(upListNow) atom = reset(defaultFirstAtom());
        else if(doGoDown) atom = reset(defaultLastAtom());
        else atom = null;
        hasNext = (atom != null);
        return atom;
    }

    /**
     * Resets the iterator so that it is ready to go through its list beginning from
     * the given atom if the initiation flag is set to INCLUDE_FIRST, otherwise it
     * begins with the atom following the given one.
     * Iterations proceed until the natural termination of the sequence.
     * If atom is not in list, hasNext is set to return false.
     * 
     * @param a  the nominal first atom returned by the iterator
     * @return   the atom that will be returned with the next subsequent call to next()
     */
    protected Atom reset(Atom a) {
        if(isNeighborIterator) atom = firstNeighbor(a);
        else atom = contains(a) ? a : null; //also ensures that a is not null
        
        //these are needed only if direction == BOTH
        setAtom = a;
        terminator2 = defaultFirstAtom();
        
        if(upListNow) terminator = defaultLastAtom();
        else terminator = doGoDown ? terminator2 : a;
        
  //      terminator = doGoDown ? defaultFirstAtom() : (upListNow ? defaultLastAtom() : a);
        hasNext = (atom != null);
        return atom;
    }

    /**
     * Resets iterator in reference to the given atoms.  Initiation of the iteration is
     * as given by a call to reset(first), and termination of the iteration is done in reference
     * to the second given atom.  If this atom is among the iterates of this iterator, iteration
     * will terminate when this atom is reached and returned by next()
     * If either first or last are not among the iterates, hasNext is set to return false.
     *
     * @param first the nominal first atom returned by the iterator
     * @param last  the nominal last atom returned by the iterator
     * @return      the atom that will be returned with the first call to next()
     */
     protected Atom reset(Atom first, Atom last) {
        if(isNeighborIterator) {
            if(contains(last)) { //won't work right is last is downlist of first
                atom = firstNeighbor(first);
                terminator = last;
            }
            else atom = null;
        }
        else {
            if(isOrdered(first,last) && upListNow) {
                doGoDown = false;//go through list only once
                reset(first);
                terminator = last;
            }
            else if(isOrdered(last, first) && doGoDown) {
                upListNow = false;
                reset(first);
                terminator = last;
            } 
            else atom = null;
        }
        hasNext = atom != null;
        return atom;
    }//end of reset(Atom, Atom)
    
}//end of AtomIteratorAbstract
    
