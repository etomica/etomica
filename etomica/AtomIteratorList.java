package etomica;

/**
 * Atom iterator that traverses the elements of an AtomList.
 * Configurable to permit iteration up and/or down list, beginning
 * with any specified atom in list, and ending at any specified Tab.
 *
 * @author David Kofke
 */
 
//TODO try to enable handling of possibility of setting first as an atom linker that is not in the list

 /* History of changes
  * 09/01/02 (DAK) modified nextLinker method to properly handle case of NEITHER direction
  * 08/23/04 (DAK, AS, KB) overhauled with revision of iterators
  */
public final class AtomIteratorList implements AtomIterator {
    
    private AtomList list;
    
	private AtomLinker first;//first atom or tab for iteration
	private AtomLinker.Tab terminator;//tab indicating end of iteration
	private AtomLinker next;//holds next atom to return on call to next()
    private boolean upList;
    private final Atom[] atoms = new Atom[1];
    
    /**
     * Constructs a new iterator using an empty list as its basis for iteration.
     */
	public AtomIteratorList() {
	    this(new AtomList());
	}
	/**
	 * Loops through the given iterator as currently reset and constructs a
	 * new list of atoms from its iterates, and constructs a new iterator
	 * using this new list as its basis for iteration.  Iterator is reset
	 * for iteration upon construction (constructor calls reset()).
	 */
	public AtomIteratorList(AtomIterator iterator) {
	    this(new AtomList(iterator));
	}
	/**
	 * Constructs a new iterator using the given list as its basis for iteration.
	 * Iterator is reset for iteration upon constrcution (constructor calls reset()).
	 */
	public AtomIteratorList(AtomList list) {
	    setList(list);
	    reset();
	}
	
	/**
	 * Makes a linked-list copy of the given list and constructs
	 * an iterator using the copy as a basis.  Useful to iterate over
	 * a list of atoms while doing operations that might change
	 * their order in the original list.
	 */
	public static AtomIteratorList makeCopylistIterator(AtomList list) {
	    return new AtomIteratorList(new AtomIteratorList(list));
	}
	
	
	/**
	 * Returns true if this iterator has another iterate, false otherwise.
	 */
	public boolean hasNext() {return next.atom != null;}
	    
    /**
     * Sets the given list of atoms as the basis for iteration.  The atoms
     * returned by this iterator will be those from the given list.  A subsequent
     * call to one of the reset methods is required before iterator is ready 
     * for iteration (until then, hasNext is false).  Sets basis to null.
     */
    public void setList(AtomList newList) {
        list = (newList != null) ? newList : new AtomList();
        next = terminator = list.header;
        first = list.header.next;
    }
    
    public AtomList getList() {
    	return list;
    }

    /**
     * Resets the iterator using the current values of first, terminator, and upList.  
     */
    public void reset() {
    	next = first;
    	while(next.atom == null && next != terminator && next != list.header) {
            next = upList ? next.next : next.previous;
            if(terminator == null) break;
        }
    }

	/**
     * Performs action on all atoms as prescribed in most recent call to reset.
     * Set of atoms for this method is same as that which would be given
     * by a hasNext/next loop.
     */
    public void allAtoms(AtomsetActive action){
    	AtomLinker.Tab header = list.header;
    	for(AtomLinker link=first; link!=terminator && link!=header; link=(upList ? link.next : link.previous)) {
    		if(link.atom != null) {
    			atoms[0] = link.atom;
    			action.actionPerformed(atoms);
    		}
    		else if(terminator == null) break;
    	}       
    }//end of allAtoms
    
    /**
     * Sets iteration to be in the given direction and unsets iterator.
     */
    public void setIterationDirection(IteratorDirective.Direction direction) {
        upList = (direction == IteratorDirective.UP);
		unset();
    }
    
    public IteratorDirective.Direction getIterationDirection() {
    	return upList ? IteratorDirective.UP : IteratorDirective.DOWN;
    }
    
    /**
     * Resets to begin with the given atom linker.  
     * Does not check that the linker is an iterate of this iterator.
     * If given linker is null, sets first to header of list.
     */
    public void setFirst(AtomLinker first) {
    	this.first = (first == null) ? list.header : first;
        unset();
    }
    
    /**
     * Resets in reference to the given atom.  Finds the atom in the list and
     * calls reset(AtomLinker) in reference to its linker.  If atom is not in list,
     * or is null, hasNext will be false.
     */
    public void setFirst(Atom atom) {
        setFirst(list.entry(atom));
    }

    public AtomLinker getFirst() {return first;}
    
    /**
     * Resets for new iteration, beginning with the atom of the first argument.
     * If first is an index, iterator is advanced to begin with the
     * next atom entry.
     */
    public void setTerminator(AtomLinker.Tab terminator) {
    	if((terminator != null) && (terminator.list != this.list)) throw new IllegalArgumentException("Error in setting terminator as an element not in the list set for iteration");
        this.terminator = terminator;
        unset();
    }
    
    public AtomLinker.Tab getTerminator() {
    	return terminator;
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
    public void set(AtomLinker first, AtomLinker.Tab terminator, IteratorDirective.Direction direction) {
        this.first = first;
        this.terminator = terminator;
        setIterationDirection(direction);
    }
    
    /**
     * Sets iterator such that hasNext() will return false.
     */
    public void unset() {
        next = list.header;
    }

    /**
     * Returns true if the given atom is in the list of iterates, false otherwise.
     */
	public boolean contains(Atom[] atom){
        if(first == list.header && terminator == list.header) return list.contains(atom[0]);
		AtomsetActiveDetect detector = new AtomsetActiveDetect(atom[0]);
		allAtoms(detector);
		return detector.detectedAtom();
	}
	
	/**
	 * Returns the total number of iterates that can be returned by this iterator, for
	 * its current list basis.
	 */
	public int size() {
		if(first == list.header && terminator == list.header) return list.size();
		AtomsetActiveCount counter = new AtomsetActiveCount();
		allAtoms(counter);
		return counter.callCount();

	}
	
	public Atom[] next() {
		atoms[0] = nextLinker().atom;
		return atoms;
	}
		
	/**
	 * Returns the next iterate.
	 */
    public Atom nextAtom() {
        return nextLinker().atom;
    }
    
    /**
     * Returns the next atom in the list without advancing the iterator.
     */
    public Atom[] peek() {
    	atoms[0] = next.atom;
        return atoms;
    }
    
    public final int nBody() {return 1;}
    
    public AtomLinker nextLinker() {
    	if(next.atom == null) return next;//prevent call to nextLinker from advancing past terminator
        AtomLinker nextLinker = next;
        next = upList ? next.next : next.previous;
        while(next.atom == null) {
            //if terminator is null we stop at the first encounter of a Tab linker
            //otherwise stop only if Tab linker is the specified terminator or the header (which could be encountered before terminator, if different
            if(terminator == null || next == terminator || next == list.header) break;//check against header also, in case it is not the terminator but it is reached first
            next = upList ? next.next : next.previous;
        }
        return nextLinker;
    }//end of nextLinker
    
    /**
     * Method to test and demonstrate use of class.
     */
    public static void main(String[] args) {
        
        Simulation sim = new Simulation();
        SpeciesSpheresMono species = new SpeciesSpheresMono();
        species.setNMolecules(10);//tested also for 0 and 1 molecule
        Phase phase = new Phase();
        sim.elementCoordinator.go();
//        AtomList atomList = phase.speciesMaster.atomList;
        AtomList atomList = ((AtomTreeNodeGroup)phase.getAgent(species).node).childList;
        
        AtomIteratorList iterator = new AtomIteratorList(atomList);
        Atom first = atomList.getFirst();
        Atom last = atomList.getLast();
        Atom middle = null;
        try {
            middle = atomList.get(atomList.size()/2);//exception thrown if list is empty
        } catch(IndexOutOfBoundsException ex) {}
        
//        iterator.setSkipFirstAtom(true);
        IteratorDirective.testSuite(iterator, first, middle, last);
    }

}//end of AtomIteratorList

