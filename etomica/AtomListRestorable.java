package etomica;

/**
 * An extension of the AtomList class the permits the list to be restored
 * to an earlier condition.  List tracks all additions/removals to it,
 * and can be instructed to undo all such changes that occurred from a 
 * particular time.  Upon restore(), list will contain all and only those atoms
 * that were in list at construction or at last call to clearMemory().
 * Restoration does not necessarily put atoms back in original order;
 * only the contents of the list is guaranteed to be as before.
 *
 * @author David Kofke
 * 02.03.21 
 */
 
 public class AtomListRestorable extends AtomList {
 
    private final AtomList memoryList = new AtomList();
    
    public AtomListRestorable() {
        super();
    }
    
    public AtomListRestorable(AtomIterator iterator) {
        super(iterator);
        clearMemory();
    }
    
    /**
     * Adds first linker before the second one, and notes its removal in the
     * memory of changes.  All add methods inherited from AtomList go through
     * this method, and all changes made with them will be recorded similarly.
     */
	public AtomLinker addBefore(AtomLinker newAtomLinker, AtomLinker e) {
	    ((AtomLinkerRestorable)newAtomLinker).wasAdded = true;
	    memoryList.add(((AtomLinkerRestorable)newAtomLinker).shadowLinker);
	    return super.addBefore(newAtomLinker, e);
	}
	
	/**
	 * Removes linker from list and notes its removal in memory of changes.
	 * All remove methods inherited from AtomList go through this method, and
	 * all changes made with them will be recorded similarly.
	 */
	public void remove(AtomLinker e) {
	    super.remove(e);
	    ((AtomLinkerRestorable)e).wasAdded = false;
	    memoryList.add(((AtomLinkerRestorable)e).shadowLinker);
	}
	
	/**
	 * Puts contents of list to condition equal to that at last call to
	 * clearMemory().  Order of atoms in list is not necessarily as before.
	 */
	public void restore() {
	    while(memoryList.size() > 0) {
	        AtomLinkerRestorable next = (AtomLinkerRestorable)memoryList.lastEntry();
	        if(next.shadowLinker.wasAdded) {//entry was added to list; remove it
	            super.remove(next.shadowLinker);
	        } else {//was removed from list; put it back
	            super.addBefore(next.shadowLinker, header);//must use this form to add; all others end up calling this subclass's version of addBefore
	        }
	        memoryList.remove(next);
	    }
	}
	
	/**
	 * Clears the memory that records added/removed atoms.
	 */
	public void clearMemory() {
	    memoryList.clear();
	}
	
	/**
	 * Clears the list but not the memory of changes, so can be undone with restore.
	 */
	public void clear() {
	    while(size() > 0) remove(lastEntry());
	 /*   super.clear();
	    clearMemory();*/
	}
    
    /**
     * Overrides the superclass method to make a linker appropriate for
     * maintenance of the memory of changes.
     */
    protected AtomLinker makeLinker(Atom atom) {
        return new AtomLinkerRestorable(atom);
    }
    
    /**
     * Extension of AtomLinker that holds a flag used to indicate
     * if atom was added or removed from list.  Also holds a "shadow"
     * linker of the same type.  Linker and its partner have
     * each other as shadow linkers; one linker is used to for this class' linked list,
     * and the other is used for maintenance of the memory of changes.
     */
    private class AtomLinkerRestorable extends AtomLinker {
        boolean wasAdded = true;
        final AtomLinkerRestorable shadowLinker;
        AtomLinkerRestorable(Atom atom) {
            super(atom);
            shadowLinker = new AtomLinkerRestorable(atom, this);
        }
        private AtomLinkerRestorable(Atom atom, AtomLinkerRestorable shadowLinker) {
            super(atom);
            this.shadowLinker = shadowLinker;
        }
    }//end of AtomLinkerRestorable
    
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
        
        list.remove(7);
        list.remove(2);
        System.out.println("Removed elements 7 and 2");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        if(pauseForInput) IteratorDirective.pauseForInput();
        
        list.add(species.factory.makeAtom());
        list.add(species.factory.makeAtom());
        list.add(species.factory.makeAtom());
        System.out.println("Added 3 atoms");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        if(pauseForInput) IteratorDirective.pauseForInput();
        
        list.restore();
        System.out.println("Restored");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        if(pauseForInput) IteratorDirective.pauseForInput();
        
        list.restore();
        System.out.println("Restored again (should be no change)");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        if(pauseForInput) IteratorDirective.pauseForInput();
        
        list.clear();
        System.out.println("Cleared list");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        if(pauseForInput) IteratorDirective.pauseForInput();
        
        list.restore();
        System.out.println("Restored");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        if(pauseForInput) IteratorDirective.pauseForInput();
 
        list.add(species.factory.makeAtom());
        list.add(species.factory.makeAtom());
        System.out.println("Added 2 atoms");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        if(pauseForInput) IteratorDirective.pauseForInput();
        
        list.clearMemory();
        list.add(species.factory.makeAtom());
        System.out.println("Cleared memory and added 1 atom");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        if(pauseForInput) IteratorDirective.pauseForInput();
        
        list.restore();
        System.out.println("Restored (should have two extra atoms)");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        if(pauseForInput) IteratorDirective.pauseForInput();
        
        list.clear();
        list.clearMemory();
        System.out.println("Cleared list and memory");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        if(pauseForInput) IteratorDirective.pauseForInput();
        
        list.restore();
        System.out.println("Restored; should still be empty");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        if(pauseForInput) IteratorDirective.pauseForInput();
        
        list.add(species.factory.makeAtom());
        list.add(species.factory.makeAtom());
        System.out.println("Added 2 atoms");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        if(pauseForInput) IteratorDirective.pauseForInput();
        
        list.restore();
        System.out.println("Restored; should be empty");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        if(pauseForInput) IteratorDirective.pauseForInput();
        
        list.restore();
        System.out.println("Restored again (should be no change)");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        if(pauseForInput) IteratorDirective.pauseForInput();
        
    }//end main
 }//end AtomListRestorable