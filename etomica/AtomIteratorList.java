package etomica;

/**
 * Atom iterator that traverses the elements of an AtomList.
 * IteratorDirective direction can set whether iterator goes up or down list.
 * Can be set or as a neighbor iterator, in which case the first indicated atom
 * is skipped.
 *
 * @author David Kofke
 */
public final class AtomIteratorList implements AtomIterator {
    
	private AtomLinker next, start, header;
    private AtomList list;
    protected boolean hasNext;
    protected boolean upListNow, doGoDown;
    private boolean isAsNeighbor;
    private IteratorDirective.Direction direction;
        
	AtomIteratorList() {
	    this(AtomList.NULL);
	}
	public AtomIteratorList(AtomList list) {
	    setBasis(list);
	}
	
	public void reset(int index) {
	    if (index < 0 || index > list.size())
		throw new IndexOutOfBoundsException("Index: "+index+
						    ", Size: "+list.size());
	    if (index < list.size()/2) {
		    next = header.next;
		    for (int nextIndex=0; nextIndex<index; nextIndex++) next = next.next;
	    } else {
		    next = header;
		    for (int nextIndex=list.size(); nextIndex>index; nextIndex--) next = next.previous;
	    }
	    hasNext = (next != header);
	}

	public boolean hasNext() {return hasNext;}
	
    public void setBasis(Atom dummyAtom){ /* no implementation */}
    public void setBasis(AtomList list) {
        this.list = list;
        if(list == null) {hasNext = false; return;}
        header = list.header;
        reset();
    }
    public Atom getBasis(){return null;}//no implementation
    
    public void setAsNeighbor(boolean b){
        isAsNeighbor = b;
    }
    
    public void allAtoms(AtomAction action){
        for (AtomLinker e = header.next; e != header; e = e.next) action.actionPerformed(e.atom);
    }
    
    /**
     * Sets iterator so that it is ready to go up its entire list of iterates.
     */
    public Atom reset() {
        isAsNeighbor = false;
        direction = IteratorDirective.UP;
        applyDirection();
        return reset(header.next);
    }
    
    /**
     * Resets using direction and atom specified in directive.
     */
    public Atom reset(IteratorDirective id){
        direction = id.direction();
        applyDirection();
        switch(id.atomCount()) {
            case 0:  return reset(upListNow ? header.next : (doGoDown ? header.previous : null)); 
            case 1:  return reset(id.atom1()); 
            default: hasNext = false; 
            return null;
        }
    }
    
    /**
     * Resets in reference to the given atom linker.  If not isAsNeighbor, the
     * given linker (or its atom) is the next iterate; if isAsNeighbor, an adjacent
     * iterate (up or down as given by the current direction) is the next iterate.
     * Does not check that the linker is an iterate of this iterator.
     */
    public Atom reset(AtomLinker linker) {
        if(linker == null) {
            hasNext = false;
            return null;
        }
        if(isAsNeighbor) {
            next = upListNow ? linker.next : (doGoDown ? linker.previous : null);
            if(next == null) {
                hasNext = false;
                return null;
            }
        }
        else next = linker;
        hasNext = true;
        start = next;
        return next.atom;
    }
    
    /**
     * Resets in reference to the given atom.  Finds the atom in the list and
     * calls reset(AtomLinker) in reference to its linker.  If atom is not in list,
     * hasNext is set to false.
     */
    public Atom reset(Atom atom) {
        return reset(list.linker(atom));
    }
    
    private void applyDirection() {
        upListNow = direction.doUp();
        doGoDown = direction.doDown();
    }

    /**
     * Returns true if the given atom is in the list of iterates, false otherwise.
     */
	public boolean contains(Atom atom){
        for (AtomLinker e=header.next; e!=header; e=e.next) if(e.atom == atom) return true;
        return false;
	}
	
	/**
	 * Returns the total number of iterates that can be returned by this iterator, for
	 * its current list basis.
	 */
	public int size(){ return list.size();}

	    
    public Atom next() {
        return nextLinker().atom;
    }
    
    public AtomLinker nextLinker() {
        AtomLinker nextLinker = next;
        next = upListNow ? next.next : next.previous;
        hasNext = (next != header);
        if(!hasNext && upListNow && doGoDown) {//done going up and now prepare to go down
            next = start.previous;
            hasNext = next != header;
            upListNow = false;
        }
        return nextLinker;
    }

}//end of AtomIteratorList

