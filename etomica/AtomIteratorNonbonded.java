package etomica;

/**
 * Loops over all the atoms in a basis which are not bonded to 
 * an atom specified in the iterator directive.
 */
 
 /* History
  * 12/06/02 (DAK) deleted line in reset(IteratorDirective) that resets basis using reference atom
  */

public class AtomIteratorNonbonded implements AtomIterator {
    
    private final AtomIterator iterator;
    
    private boolean hasNext;
    protected IteratorDirective.Direction direction;
    private Atom nonbondAtom;
    private Atom nextAtom;
    
    public AtomIteratorNonbonded(Simulation sim) {
        iterator = sim.iteratorFactory.makeIntragroupNbrIterator();
    }
    
	public void all(AtomSet basis, IteratorDirective id, final AtomsetActive action) {
		 if(!(basis instanceof Atom && action instanceof AtomActive)) return;
		 all((Atom)basis, id, (AtomActive)action);
	}
    
	/**
	 * this throws because it can't work!!!!
	 */
	public void all(Atom basis, IteratorDirective id, final AtomActive action) {
		if(basis == null || basis.node.isLeaf() || action == null) return;
//		iterator.all(basis, id, action);
		throw new RuntimeException("Method all not implemented in AtomIteratorNonbonded");
	}

    public boolean hasNext() {return hasNext;}

    /**
     * Puts iterator in a state in which hasNext is false.
     */
    public void unset() {hasNext = false;}

    public boolean contains(Atom atom) {
        return(iterator.contains(atom) && !Bond.areBonded(nonbondAtom, atom));
    }
    
    public void reset(IteratorDirective id) {
        nonbondAtom = id.atom1();
// commented 12-06-02        if(nonbondAtom != null) setBasis(nonbondAtom.node.parentGroup());
        iterator.reset(id);
        next();
    }
    
    public void reset() {
        iterator.reset();
        next();
    }
    
    public Atom next() {
        Atom next = nextAtom;
        hasNext = false; nextAtom = null;
        while(iterator.hasNext()) {//loop until another nonbonded atom is found, or iterator expires
            nextAtom = iterator.next();
            if(!Bond.areBonded(nextAtom, nonbondAtom)) {
                hasNext = true;
                break;
            }
        }
        return next;
    }
    
    //not implemented
    public void allAtoms(AtomAction act) {}
    
    public void setBasis(Atom atom) {
        iterator.setBasis(atom);
    }
    
    public Atom getBasis() {return iterator.getBasis();}
    
    /**
     * Do not use: not correctly implemented.
     */
    public int size() {
        throw new RuntimeException("AtomIteratorNonbonded.size() not implemented");
    }   
}//end of AtomIteratorNonbonded