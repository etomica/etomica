package etomica;

/**
 * Set of AtomIterators that are selected according to the direction indicated
 * by the iteratorDirective.
 */
public final class AtomIteratorBundle extends AtomIterator {
        
    AtomIterator current, up, down, singlet;
    Phase phase;
    private boolean neighborIterator;
        
    public AtomIteratorBundle(Phase p) {
        phase = p;
        makeAtomIterators();
        
        //register to be informed if iterator factory in phase is changed
        phase.iteratorFactoryMonitor.addObserver(new java.util.Observer() {
	        public void update(java.util.Observable o, Object arg) {
	            makeAtomIterators();
	        }
	    });
	}//end of constructor
	    
    private void makeAtomIterators() {
        up = phase.iteratorFactory().makeAtomIteratorUp();
        down = phase.iteratorFactory().makeAtomIteratorDown();
        singlet = new AtomIterator.Singlet();
        current = singlet;
    }//end of makeAtomIterators
        
    public void reset(IteratorDirective id) {
        IteratorDirective.Direction direction = id.direction();
            
        if(direction == IteratorDirective.UP) current = up;
        else if(direction == IteratorDirective.SINGLET) current = singlet;
        else current = down;
            
        current.reset(id);
    }
    public Atom next() {return current.next();}
        
    public Atom reset(Atom a) {return current.reset(a);}

    /**
    * Resets the iterator, so that it is ready to go through its list again.
    */
    public Atom reset() {return current.reset();}

    /**
    * Resets iterator so that it loops through the given atoms, inclusive.
    */
    public Atom reset(Atom first, Atom last) {return current.reset(first, last);}
    /**
    * Performs the given Action on each atom in the list in sequence.
    * 
    * @param act
    * @see Atom.Action
    */
    public void allAtoms(AtomAction act) {current.allAtoms(act);}
}//end of Bundle
            
