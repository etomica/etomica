package etomica;

/**
 * Basic interface for iterating over pairs of atoms.
 *
 * @author David Kofke
 */
public abstract class AtomPairIterator implements AtomSetIterator, java.io.Serializable {
	
	public void all(AtomSet basis, IteratorDirective id, final AtomSetActive action) {
	 if(basis == null || !(action instanceof AtomPairActive)) return;
	 switch(basis.atomCount()) {
		case 1: all((Atom)basis, id, (AtomPairActive)action); break;
		case 2: all((AtomPair)basis, id, (AtomPairActive)action); break;
	 }
	}

	public abstract void all(Atom basis, IteratorDirective id, AtomPairActive action);
	
	public abstract void all(AtomPair basis, IteratorDirective id, AtomPairActive action);
    
    public abstract void setBasis(Atom a1, Atom a2);
    
    /**
     * Returns the number of pairs that would be given by this iterator
     * after a call to the no-argument reset() method.
     */
    public abstract int size();        
    
    public abstract boolean hasNext();
    
    public abstract void reset(IteratorDirective id);
    
    
   /**
     * Resets the iterator, so that it is ready to go through all of its pairs.
     */
    public abstract void reset();
        
        
    public abstract AtomPair next();

    /**
     * Performs the given action on all pairs returned by this iterator.
     */
    public abstract void allPairs(AtomPairAction act);
    
    public static final AtomPairIterator NULL = new Null();
    static final class Null extends AtomPairIterator {
        public void setBasis(Atom a1, Atom a2) {}
        public int size() {return 0;}                
        public boolean hasNext() {return false;}       
        public void reset(IteratorDirective id) {}
        public void reset() {}
        public void reset(Atom atom) {}
        public AtomPair next() {return null;}
        public void allPairs(AtomPairAction act) {}
		public void all(Atom basis, IteratorDirective id, AtomPairActive act) {}
        public void all(AtomPair basis, IteratorDirective id, AtomPairActive act) {}
    }//end of Null
    
}  //end of class AtomPairIterator
    
