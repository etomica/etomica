package etomica;

/**
 * Set of AtomIterators that are selected according to the direction indicated
 * by the iteratorDirective.
 */
public final class AtomIteratorBundle implements AtomIterator {
        
    AtomIterator current, up, down, singlet, both;
    Phase phase;
    private AtomIterator.Initiation initiation;
    private boolean neighborIterator;
        
    public AtomIteratorBundle(Phase p) {
        this(p, AtomIterator.INCLUDE_FIRST);
    }
    public AtomIteratorBundle(Phase p, AtomIterator.Initiation init) {
        phase = p;
        initiation = init;
        makeDefaultIterators();
        
        //register to be informed if iterator factory in phase is changed
        phase.iteratorFactoryMonitor.addObserver(new java.util.Observer() {
	        public void update(java.util.Observable o, Object arg) {
	            makeDefaultIterators();
	        }
	    });
	}//end of constructor
	    
    private void makeDefaultIterators() {
        up = phase.iteratorFactory().makeAtomIteratorUp(initiation);
        down = phase.iteratorFactory().makeAtomIteratorDown(initiation);
        both = phase.iteratorFactory().makeAtomIteratorUpDown(initiation);
        singlet = new AtomIteratorSinglet(initiation);
        current = singlet;
    }//end of makeAtomIterators
    
    public void setIterators(AtomIterator up, AtomIterator down,
                             AtomIterator singlet, AtomIterator both) {
        if(up != null) this.up = up;
        if(down != null) this.down = down;
        if(singlet != null) this.singlet = singlet;
        if(both != null) this.both = both;
    }
    
    public boolean hasNext() {return current.hasNext();}
        
    public void reset(IteratorDirective id) {
        IteratorDirective.Direction direction = id.direction();
            
        if(direction == IteratorDirective.UP) current = up;
        else if(direction == IteratorDirective.NEITHER) current = singlet;
        else if(direction == IteratorDirective.DOWN) current = down;
        else current = both;
            
        current.reset(id);
    }
    
    public boolean contains(Atom atom) {return current.contains(atom);}
    
    public Atom next() {return current.next();}
        
    /**
    * Performs the given Action on each atom in the list in sequence.
    * 
    * @param act
    * @see Atom.Action
    */
    public void allAtoms(AtomAction act) {current.allAtoms(act);}
}//end of AtomIteratorBundle
            
