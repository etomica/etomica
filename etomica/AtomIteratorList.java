package etomica;

/**
 * Atom iterator that traverses the elements of an AtomList.
 * IteratorDirective direction can set to cause iterator to go up or down list.
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
	
	/**
	 * Sets the childList of the given atom as the basis for iteration.
	 * If atom is a leaf, hasNext is false.
	 */
    public void setBasis(Atom atom){
        if(atom.node.isLeaf()) unset();
        else setBasis(((AtomTreeNodeGroup)atom.node).childList);
    }
    public void setBasis(AtomList list) {
        this.list = list;
        if(list == null) header = nullLinker;
        else header = list.header;
        next = terminator = header;
    }
    public Atom getBasis(){return null;}//no implementation

    /**
     * Performs action on all atoms.
     */
     //needs development to perform action according to settings of terminator etc.
    public void allAtoms(AtomAction action){
        for (AtomLinker e = header.next; e != header; e = e.next) 
            if(e.atom != null) action.actionPerformed(e.atom);
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
    
    /**
     * Resets to begin iteration in the given direction, stopping when the specified tab
     * is encountered.  If direction is UP, iteration proceeds up list and ends when terminator
     * or header is encountered; likewise if direction is DOWN.  If direction is BOTH, 
     * proceeds up list from starting point until encountering header or terminator, 
     * and then down it from the starting point, again until encountering header or terminator.
     * If terminator is null, iteration halts in each direction when any tab (or the header) 
     * is encountered.<br>
     * To iterate completely in either or both directions (ignoring all tabs), use the 
     * reset(AtomLinker, Direction) method.<br>
     * Up iteration always begins by returning the given first linker
     * (unless it is a tab); down iteration always begins with the linker before the given
     * first one.
     */
    public Atom reset(AtomLinker first, AtomLinker.Tab terminator, IteratorDirective.Direction direction) {
        if(first == header || first == null) {
            next = header;
            return null;
        }
        applyDirection(direction);
        this.terminator = terminator;
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
    
    public AtomLinker nextLinker() {
        AtomLinker nextLinker = next;
        next = upListNow ? next.next : next.previous;
        while(next.atom == null) {
            //if terminator is null we stop at the first encounter of a Tab linker
            //otherwise stop only if Tab linker is the specified terminator
            if(terminator == null || next == header || next == terminator) {
                if(upListNow && doGoDown) {//done going up and now prepare to go down
                    next = start.previous;
                    upListNow = false;
                } else {
                    break;
                }
            }
            else next = upListNow ? next.next : next.previous;
        }
        return nextLinker;
    }//end of nextLinker

}//end of AtomIteratorList

