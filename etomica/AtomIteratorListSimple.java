package etomica;

/**
 * Lightweight version of AtomIteratorList.  Iterates uplist only, from
 * beginning to end of list.  Iterator functions correctly in situations where
 * elements are removed from list after they are returned by iterator.
 *
 * @author David Kofke
 */

/* History
 * 02/21/03 (DAK) added all(AtomList...) method
 * 
 */
public final class AtomIteratorListSimple implements AtomIterator {
    
    private AtomList list;
	private AtomLinker next;
        
	public AtomIteratorListSimple() {
	    this(AtomList.NULL);
	}
	public AtomIteratorListSimple(AtomList list) {
	    setBasis(list);
	    reset();
	}
    
	public void all(Atom basis, IteratorDirective dummy, final AtomActive action) {
		if(basis==null || basis.node.isLeaf()) return;
		all(((AtomTreeNodeGroup)basis.node).childList, dummy, action);
	}
	
	public void all(AtomList basisList, IteratorDirective dummy, final AtomActive action) {
		if(basisList == null || action == null) return;
		final AtomLinker header = basisList.header;
		for(AtomLinker e=header.next; e!=header; e=e.next) {
			if(e.atom != null) action.actionPerformed(e.atom);
		}
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
        if(list == null) list = AtomList.NULL;
        this.list = list;
        next = list.header;
    }
    public Atom getBasis(){return null;}//no implementation
 
     /**
     * Performs action on all atoms.
     */
    public void allAtoms(AtomAction action){
        for (AtomLinker e = list.header.next; e != list.header; e = e.next) 
            if(e.atom != null) action.actionPerformed(e.atom);
    }
    
    /**
     * Sets iterator so that it is ready to go up its entire list of iterates.
     */
    public Atom reset() {
        next = list.header.next;
        while(next.atom == null && next != list.header) next = next.next;
        return next.atom;
    }
    
    /**
     * Resets, ignoring given directive (same as reset()).
     */
    public Atom reset(IteratorDirective id) {return reset();}
    
    /**
     * Sets iterator such that hasNext() will return false.
     */
    public void unset() {next = list.header;}
    
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
	public int size() {return list.size();}

	    
    public Atom next() {
        return nextLinker().atom;
    }
    public AtomSet nextSet() {
    	return next();
    }
    
    /**
     * Returns the next atom in the list without advancing the iterator.
     */
    public Atom peek() {
        return next.atom;
    }
    
    private AtomLinker nextLinker() {
        AtomLinker nextLinker = next;
        next = next.next;
        while(next.atom == null && next != list.header) {next = next.next;}
        return nextLinker;
    }//end of nextLinker
    
	/**
	 * Invokes all(Atom, IteratorDirective, AtomActive) method of this
	 * class, using given arguments if they are instances of the appropriate
	 * classes. Otherwise returns without throwing any exception.
	 * @see etomica.AtomSetIterator#all(AtomSet, IteratorDirective, AtomSetActive)
	 */
	public void all(AtomSet basis, IteratorDirective id, final AtomSetActive action) {
		 if(!(basis instanceof Atom && action instanceof AtomActive)) return;
		 all((Atom)basis, id, (AtomActive)action);
	}
    
    public static void main(String[] args) {
        Simulation sim = new Simulation();
        Phase phase = new Phase();
        SpeciesSpheresMono species = new SpeciesSpheresMono();
        species.setNMolecules(10);
        sim.elementCoordinator.go();
        
        boolean pauseForInput = true;
        
        AtomListRestorable list = new AtomListRestorable(phase.makeMoleculeIterator());
        AtomIteratorListSimple iterator = new AtomIteratorListSimple(list);
        
        System.out.println("Original list");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        if(pauseForInput) IteratorDirective.pauseForInput();
        
        System.out.println("Removing each element from list as iterated");
        iterator.reset();
        while(iterator.hasNext()) {
            Atom atom = iterator.next();
            System.out.println(atom.toString());
            list.remove(atom);
        }
        if(pauseForInput) IteratorDirective.pauseForInput();
        
        System.out.println("Empty list");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        if(pauseForInput) IteratorDirective.pauseForInput();
    }//end main

}//end of AtomIteratorListSimple

