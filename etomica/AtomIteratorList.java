package etomica;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.action.AtomsetDetect;

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
public final class AtomIteratorList implements AtomIteratorListDependent, AtomsetIteratorDirectable {
    
    private AtomList list;
    
	private AtomLinker first;//first atom or tab for iteration
	private int terminatorType;//bitmask to match against for terminators
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
	 * Loops through the given iterator as currently conditioned, constructs a
	 * new list of atoms from its iterates, and constructs a new iterator
	 * using this new list as its basis for iteration.  Given iterator need not be
	 * reset as passed to this constructor; iterator made by this constructor
	 * must be reset before use.
	 */
	public AtomIteratorList(AtomIterator iterator) {
	    this(new AtomList(iterator));
	}
	/**
	 * Constructs a new iterator using the given list as its basis for iteration.
	 * Iterator is conditioned to start from the first element of list and iterate
	 * up. Iterator must be reset before use.
	 */
	public AtomIteratorList(AtomList list) {
	    setList(list);
	    upList = true;
	}
	
	/**
	 * Makes a linked-list copy of the given list and constructs
	 * an iterator using the copy as a basis.  Useful to iterate over
	 * a list of atoms while doing operations that might change
	 * their order in the original list.  Iterator returned by this
	 * method must be reset before use.
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
     * returned by this iterator will be those from the given list.  After 
     * calling this method, iterator must be reset before use.  This method 
     * has no effect on iteration direction, but does set first and terminator 
     * iteration elements to be the header of the given list. If given
     * a null argument, an empty list is created for iteration. 
     */
    public void setList(AtomList newList) {
        list = (newList != null) ? newList : new AtomList();
        next = list.header;
        terminatorType = AtomLinker.Tab.HEADER_TAB;
        first = list.header;
    }
    
    
    /**
     * @return the list defining the atoms given by this iterator. 
     */
    public AtomList getList() {
    	return list;
    }

    /**
     * Set the iterator to begin iteration in its current condition (first, terminator
     * and direction most recently set).  
     */
    public void reset() {
    	next = first;
    	if(next.atom == null) next = upList ? next.next : next.previous;
        while(next.atom == null && (((AtomLinker.Tab)next).type & terminatorType) == 0) {
            next = upList ? next.next : next.previous;
        }
    }

	/**
     * Performs action on all atoms according to the current condition of
     * the iterator (first, terminator, direction). Does not require reset
     * before use.
     */
    public void allAtoms(AtomsetAction action){
    	for(AtomLinker link = (first.atom==null) ? (upList ? first.next : first.previous) : first; ; 
    			link=(upList ? link.next : link.previous)) {
    		if(link.atom != null) {
    			atoms[0] = link.atom;
    			action.actionPerformed(atoms);
    		}
    		else if((((AtomLinker.Tab)link).type & terminatorType) > 0) break;
    	}       
    }//end of allAtoms
    
    /**
     * Sets iteration to be in the given direction and unsets iterator.
     */
    public void setDirection(IteratorDirective.Direction direction) {
        upList = (direction == IteratorDirective.UP) || (direction == null);
		unset();
    }
    
    /**
     * @return the current setting for the iteration direction.
     */
    public IteratorDirective.Direction getDirection() {
    	return upList ? IteratorDirective.UP : IteratorDirective.DOWN;
    }
    
    /**
     * Resets to begin with the given atom linker and unsets iterator.  
     * Does not check that the linker is an iterate of this iterator, and
     * iterator will not function correctly if linker is not in the current list.
     * If given linker is null, sets first to header of list.
     */
    //TODO revisit the design that permits specifying linker not in list
    public void setFirst(AtomLinker first) {
    	this.first = (first == null) ? list.header : first;
        unset();
    }
    
    /**
     * Sets iteration to begin with the given atom. Finds the atom in the list and
     * calls reset(AtomLinker) in reference to its linker.  If atom is null, sets
     * iteration to begin with first element of list. If atom is not in list, behavior
     * is equivalent to setting with null atom.
     * In all cases, reset is required before using iterator.
     */
    public void setFirst(Atom atom) {
    	setFirst(list.entry(atom));//entry returns null if atom is null or is not in list
    }

    /**
     * @return the linker of the atom corresponding to the most recent call to setFirst,
     * or the list header if a first atom has not be specified.
     */
    public AtomLinker getFirst() {return first;}
    
    /**
     * Sets the tab type of the terminator.  Iteration stops when the header
     * or a tab of this type is encountered.  
     */
    public void setTerminatorType(int terminatorType) {
    	this.terminatorType = terminatorType | AtomLinker.Tab.HEADER_TAB;
    }
    
    /**
     * @return the current tab type that indicates termination of iteration
     */
    public int getTerminatorType() {
    	return terminatorType;
    }
    
    /**
     * Sets iterator such that hasNext() will return false.
     */
    public void unset() {
        next = list.header;
    }

    /**
     * Returns true if the first atom of the given array is in the list of iterates 
     * that would be given by iterator as currently conditioned; returns false otherwise.  
     * Ignores any atoms other than the first in the array.  Returns false if array does 
     * not specify an atom.  Does not require iterator be reset.
     */
	public boolean contains(Atom[] atom){
		if(atom == null || atom.length == 0) return false;
		// use AtomList's contains() if first is the header and the terminator type is header (only).
        if(first == list.header && terminatorType == AtomLinker.Tab.HEADER_TAB)	return list.contains(atom[0]);//returns false also if atom[0] is null
		AtomsetDetect detector = new AtomsetDetect(atom[0]);
		allAtoms(detector);
		return detector.detectedAtom();
	}
	
	/**
	 * Returns the total number of iterates that can be returned by this iterator, for
	 * its current list basis, and as it is currently conditioned (as given by most recently
	 * specified first and terminatorType). Does not require the iterator be reset.
	 */
	public int size() {
		// use AtomList's size() if first is the header and the terminator type is header (only).
		if(first == list.header && terminatorType == AtomLinker.Tab.HEADER_TAB) return list.size();
		AtomsetCount counter = new AtomsetCount();
		allAtoms(counter);
		return counter.callCount();
	}
	
	/**
	 * Returns the next atom from the iterator, as the first element of the returned array.
	 * Implementation of AtomsetIterator interface.
	 */
	public Atom[] next() {
		atoms[0] = nextLinker().atom;
		return atoms;
	}
		
	/**
	 * Returns the next iterate. Implementation of the AtomIterator interface.
	 */
    public Atom nextAtom() {
        return nextLinker().atom;
    }
    
    /**
     * Returns the next atom in the list without advancing the iterator.  Atom 
     * is given as the first element of the returned array.
     */
    public Atom[] peek() {
    	atoms[0] = next.atom;
        return atoms;
    }
    
    /**
     * Returns 1, indicating that the next() method returns an array with one element.
     */
    public final int nBody() {return 1;}
    
    /**
     * Returns the linker of the next atom given by the iterator.  Advances
     * iterator in doing so.
     */
    public AtomLinker nextLinker() {
    	if(next.atom == null) return next;//prevent call to nextLinker from advancing past terminator
        AtomLinker nextLinker = next;
        next = upList ? next.next : next.previous;
        // loop through any tabs that do not match the terminatorType bitmask
        while(next.atom == null && (((AtomLinker.Tab)next).type & terminatorType) == 0) {
            next = upList ? next.next : next.previous;
        }
        return nextLinker;
    }//end of nextLinker
    
    /**
     * Method to test and demonstrate use of class.
     */
    public static void main(String[] args) {
        
        Simulation sim = new Simulation();
        SpeciesSpheresMono species = new SpeciesSpheresMono(sim);
        species.setNMolecules(10);//tested also for 0 and 1 molecule
        Phase phase = new Phase(sim.space);
//        sim.elementCoordinator.go();
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
//        IteratorDirective.testSuite(iterator, first, middle, last);
    }

}//end of AtomIteratorList

