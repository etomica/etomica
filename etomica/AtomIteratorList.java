package etomica;

/**
 * Atom iterator that traverses the elements of an AtomList.
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
	    this.list = list;
	    header = list.header;
        reset();
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
        reset();
    }
    public Atom getBasis(){return null;}//no implementation
    
    public void setAsNeighbor(boolean b){
        isAsNeighbor = b;
    }
    
    public void allAtoms(AtomAction action){
        for (AtomLinker e = header.next; e != header; e = e.next) action.actionPerformed(e.atom);
    }
    
    public Atom reset() {
        isAsNeighbor = false;
        direction = IteratorDirective.UP;
        applyDirection();
        return reset(header.next);
    }
    
    public Atom reset(AtomLinker linker) {
        if(linker == null) {
            hasNext = false;
            return null;
        }
        if(isAsNeighbor) {
            next = upListNow ? linker.next : linker.previous;
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
    
    public Atom reset(IteratorDirective id){
        direction = id.direction();
        applyDirection();
        switch(id.atomCount()) {
            case 0:  return reset(list.header); 
            case 1:  return reset(id.atom1()); 
            default: hasNext = false; 
            return null;
        }
    }
    
    public Atom reset(Atom atom) {
        for (AtomLinker e=header.next; e!=header; e=e.next) {
            if(e.atom == atom) return reset(e);
        }
        hasNext = false;
        return null;
    }
    
    private void applyDirection() {
        upListNow = direction.doUp();
        doGoDown = direction.doDown();
    }

	public boolean contains(Atom atom){
        for (AtomLinker e=header.next; e!=header; e=e.next) if(e.atom == atom) return true;
        return false;
	}
	
	public int size(){ return list.size();}//

	    
    public Atom next() {
        return nextLinker().atom;
    }
    
    public AtomLinker nextLinker() {
        AtomLinker nextLinker = next;
        next = upListNow ? next.next : next.previous;
        hasNext = (next != list.header);
        if(!hasNext && upListNow && doGoDown) {//done going up and now prepare to go down
            next = start.previous;
            hasNext = next != list.header;
            upListNow = false;
        }
        return nextLinker;
    }

}//end of AtomIteratorList

