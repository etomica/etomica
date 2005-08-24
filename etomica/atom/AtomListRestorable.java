package etomica.atom;

import etomica.Phase;
import etomica.Simulation;
import etomica.SpeciesSpheresMono;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.atom.iterator.AtomIteratorListSimple;

/**
 * An extension of the AtomList class that permits the list to be restored
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
	    AtomLinkerRestorable link = (AtomLinkerRestorable)newAtomLinker;
	    link.changeCount++;
	    if(link.changeCount != 0) memoryList.add(link.shadowLinker);//maybe should be if(link.changeCount == +1)
	    else memoryList.remove(link.shadowLinker);//else if(link.changeCount == 0) no net change; remove from memory
	    return super.addBefore(newAtomLinker, e);
	}
	
	/**
	 * Removes linker from list and notes its removal in memory of changes.
	 * All remove methods inherited from AtomList go through this method, and
	 * all changes made with them will be recorded similarly.
	 */
	public void remove(AtomLinker e) {
	    AtomLinkerRestorable link = (AtomLinkerRestorable)e;
	    link.changeCount--;
	    if(link.changeCount != 0) memoryList.add(link.shadowLinker);//maybe should be if(link.changeCount == -1)
	    else memoryList.remove(link.shadowLinker);//else if(link.changeCount == 0) no net change; remove from memory
	    super.remove(e);
	}
	
	/**
	 * Puts contents of list to condition equal to that at last call to
	 * clearMemory().  Order of atoms in list is not necessarily as before.
	 */
	public void restore() {
	    while(memoryList.size() > 0) {
	        AtomLinkerRestorable next = (AtomLinkerRestorable)memoryList.lastEntry();
	        if(next.shadowLinker.changeCount > 0) {//entry was added to list; remove it
	            super.remove(next.shadowLinker);
	        } else if(next.shadowLinker.changeCount < 0) {//was removed from list; put it back
	            super.addBefore(next.shadowLinker, header);//must use this form to add; all others end up calling this subclass's version of addBefore
	        }
	        next.shadowLinker.changeCount = 0;
	        memoryList.remove(next);
	    }
	}
	
	/**
	 * Clears the memory that records added/removed atoms.
	 */
	public void clearMemory() {
	    while(memoryList.size() > 0) {
	        AtomLinkerRestorable next = (AtomLinkerRestorable)memoryList.lastEntry();
	        next.shadowLinker.changeCount = 0;
	        memoryList.remove(next);
	    }
//	    memoryList.clear();
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
    private static class AtomLinkerRestorable extends AtomLinker {
        int changeCount = 0;
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
    
    /**
     * main method to demonstrate and test use of this class.
     */
    public static void main(String[] args) {
        Simulation sim = new Simulation();
        Phase phase = new Phase(sim);
        SpeciesSpheresMono species = new SpeciesSpheresMono(sim);
        species.setNMolecules(10);
//        sim.elementCoordinator.go();
        
        boolean pauseForInput = true;
        
        AtomListRestorable list = new AtomListRestorable(new AtomIteratorAllMolecules(phase));
        AtomIteratorListSimple iterator = new AtomIteratorListSimple(list);
        
        System.out.println("Original list");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.nextAtom().toString());
        if(pauseForInput) pauseForInput();
        
        list.remove(7);
        list.remove(2);
        System.out.println("Removed elements 7 and 2");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.nextAtom().toString());
        if(pauseForInput) pauseForInput();
        
        list.add(species.moleculeFactory().makeAtom());
        list.add(species.moleculeFactory().makeAtom());
        list.add(species.moleculeFactory().makeAtom());
        System.out.println("Added 3 atoms");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.nextAtom().toString());
        if(pauseForInput) pauseForInput();
        
        list.restore();
        System.out.println("Restored");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.nextAtom().toString());
        if(pauseForInput) pauseForInput();
        
        list.restore();
        System.out.println("Restored again (should be no change)");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.nextAtom().toString());
        if(pauseForInput) pauseForInput();
        
        list.clear();
        System.out.println("Cleared list");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.nextAtom().toString());
        if(pauseForInput) pauseForInput();
        
        list.restore();
        System.out.println("Restored");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.nextAtom().toString());
        if(pauseForInput) pauseForInput();
 
        list.add(species.moleculeFactory().makeAtom());
        list.add(species.moleculeFactory().makeAtom());
        System.out.println("Added 2 atoms");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.nextAtom().toString());
        if(pauseForInput) pauseForInput();
        
        list.clearMemory();
        list.add(species.moleculeFactory().makeAtom());
        System.out.println("Cleared memory and added 1 atom");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.nextAtom().toString());
        if(pauseForInput) pauseForInput();
        
        list.restore();
        System.out.println("Restored (should have two extra atoms)");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.nextAtom().toString());
        if(pauseForInput) pauseForInput();
        
        Atom first = list.removeFirst();
        list.addFirst(species.moleculeFactory().makeAtom());
        System.out.println("Removed first and last and added a new first");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.nextAtom().toString());
        if(pauseForInput) pauseForInput();
        
        list.removeFirst();
        list.add(first);
        System.out.println("Removed new first and added original first at end");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.nextAtom().toString());
        if(pauseForInput) pauseForInput();
        
        list.restore();
        System.out.println("Restored");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.nextAtom().toString());
        if(pauseForInput) pauseForInput();
        
        list.clear();
        list.clearMemory();
        System.out.println("Cleared list and memory");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.nextAtom().toString());
        if(pauseForInput) pauseForInput();
        
        list.restore();
        System.out.println("Restored; should still be empty");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.next().toString());
        if(pauseForInput) pauseForInput();
        
        list.add(species.moleculeFactory().makeAtom());
        list.add(species.moleculeFactory().makeAtom());
        System.out.println("Added 2 atoms");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.nextAtom().toString());
        if(pauseForInput) pauseForInput();
        
        list.restore();
        System.out.println("Restored; should be empty");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.nextAtom().toString());
        if(pauseForInput) pauseForInput();
        
        list.restore();
        System.out.println("Restored again (should be no change)");
        iterator.reset();
        while(iterator.hasNext()) System.out.println(iterator.nextAtom().toString());
        if(pauseForInput) pauseForInput();
        
    }//end main

    /**
     * Halts program activity until a return is entered from the console.
     * Support method for testSuite method.
     */
    public static void pauseForInput() {
        System.out.println("Hit return to continue");
        try {
            System.in.read();
            System.in.read();
        } catch(Exception e) {}
    }

    
 }//end AtomListRestorable
