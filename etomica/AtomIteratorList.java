package etomica;

/**
 * Atom iterator that traverses the elements of an AtomList.
 * IteratorDirective direction can set to cause iterator to go up or down list.
 * Can be set or as a neighbor iterator, in which case the first indicated atom
 * is skipped.
 *
 * @author David Kofke
 */
public final class AtomIteratorList implements AtomIterator {
    
	private AtomLinker.Tab header, terminator;
    private AtomList list;
    protected boolean upListNow, doGoDown;
    private final AtomLinker.Tab nullLinker = new AtomLinker.Tab();
	private AtomLinker next = nullLinker;
	private AtomLinker start;
        
	public AtomIteratorList() {
	    this(AtomList.NULL);
	}
	public AtomIteratorList(AtomIterator iterator) {
	    this(new AtomList(iterator));
	}
	public AtomIteratorList(AtomList list) {
	    setBasis(list);
	}
	
	public void reset(int index) {
	    next = list.entry(index);
	}
	
	public void setAsNeighbor(boolean b) {
	    throw new RuntimeException("method AtomIteratorList.setAsNeighbor not implemented");
	}

	public boolean hasNext() {return next.atom != null;}
	
    public void setBasis(Atom dummyAtom){ /* no implementation */}
    public void setBasis(AtomList list) {
        this.list = list;
        if(list == null) header = nullLinker;
        else header = list.header;
        next = terminator = header;
    }
    public Atom getBasis(){return null;}//no implementation

    public void allAtoms(AtomAction action){
        for (AtomLinker e = header.next; e != header; e = e.next) action.actionPerformed(e.atom);
    }
    
    /**
     * Sets iterator so that it is ready to go up its entire list of iterates.
     */
    public Atom reset() {
        return reset(header.next, header, IteratorDirective.UP);
    }
    
    /**
     * Resets using direction and atom specified in directive.
     */
    public Atom reset(IteratorDirective id){
        applyDirection(id.direction());
        switch(id.atomCount()) {
            case 0:  return reset(upListNow ? header.next : (doGoDown ? header.previous : null)); 
            case 1:  return reset(id.atom1()); 
            default: next = header; 
            return null;
        }
    }
    
    /**
     * Resets in reference to the given atom linker.  
     * Does not check that the linker is an iterate of this iterator.
     */
    public Atom reset(AtomLinker first) {
        return reset(first, header, IteratorDirective.UP);
    }
    
    public Atom reset(AtomLinker first, IteratorDirective.Direction direction) {
        return reset(first, header, direction);
    }
    
    /**
     * Resets for new iteration, beginning with the atom of the first argument.
     * If first is an index, iterator is advanced to begin with the
     * next atom entry.
     */
    public Atom reset(AtomLinker first, AtomLinker.Tab last) {
        return reset(first, last, IteratorDirective.UP);
    }
    
    //if last == null, iterate until first Tab is encountered
    public Atom reset(AtomLinker first, AtomLinker.Tab last, IteratorDirective.Direction direction) {
        if(first == header || first == null) {
            next = header;
            return null;
        }
        applyDirection(direction);
        terminator = last;
        next = start = first;
        if(next.atom == null) nextLinker();//nextLinker keeps iterating until entry with atom is found
        return next.atom;
    }
    
    /**
     * Resets in reference to the given atom.  Finds the atom in the list and
     * calls reset(AtomLinker) in reference to its linker.  If atom is not in list,
     * hasNext will be false.
     */
    public Atom reset(Atom atom) {
        return reset(list.entry(atom), header);
    }
    
    /**
     * Sets iterator such that hasNext() will return false.
     */
    public void unset() {next = header;}
    
    private void applyDirection(IteratorDirective.Direction direction) {
        upListNow = direction.doUp();
        doGoDown = direction.doDown();
    }

    /**
     * Returns true if the given atom is in the list of iterates, false otherwise.
     */
	public boolean contains(Atom atom){
        return list.contains(atom);
	}
	
	/**
	 * Returns the total number of iterates that can be returned by this iterator, for
	 * its current list basis.
	 */
	public int size(){ return list.size();}

	    
    public Atom next() {
        return nextLinker().atom;
    }
    
    /**
     * Returns the next atom in the list without advancing the iterator.
     */
    public Atom peek() {
        return next.atom;
    }
    
/*    public AtomLinker nextLinker() {
        AtomLinker nextLinker = next;
        next = upListNow ? next.next : next.previous;
        hasNext = (next != header);
        if(!hasNext && upListNow && doGoDown) {//done going up and now prepare to go down
            next = start.previous;
            hasNext = next != header;
            upListNow = false;
        }
        return nextLinker;
    }*/

    public AtomLinker nextLinker() {
        AtomLinker nextLinker = next;
        next = upListNow ? next.next : next.previous;
        while(next.atom == null) {
            //if terminator is null we stop at the first encounter of a Tab linker
            //otherwise stop only if Tab linker is the specified terminator
            if(terminator == null || next == terminator) break;
            else next = upListNow ? next.next : next.previous;
        }
        //need to decide if we're going to permit iterator to do both directions, or just up/down
/*        hasNext = (next != header);
        if(!hasNext && upListNow && doGoDown) {//done going up and now prepare to go down
            next = start.previous;
            hasNext = next != header;
            upListNow = false;
        }*/
        return nextLinker;
    }

}//end of AtomIteratorList

