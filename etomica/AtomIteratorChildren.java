package etomica;

/**
* Iterates over all the children of a given atom group.
*
* @author David Kofke
*/

public final class AtomIteratorChildren implements AtomIterator {
    
    /**
     * Indicates if another iterate is forthcoming.
     */
    public boolean hasNext() {return next != null;}
    
    /**
     * True if the parent group of the given atom is the current basis for the iterator.
     * False otherwise, or if atom or basis is null.
     */
    public boolean contains(Atom atom) {
        return atom != null && basis != null && atom.node.parentGroup() == basis;
    }
    
    /**
     * Does reset if relation (preceeds) of atom in iterator directive and basis of this
     * iterator are compatible with the direction in iterator directive.  Otherwise
     * resets to hasNext false.
     */
    public Atom reset(IteratorDirective id) {
        if(id.direction() == IteratorDirective.BOTH 
            || id.atomCount() == 0 
            || id.atom1() == basis) return reset();
        else {
            boolean before = id.atom1().seq.preceeds(basis);
            if( (before && id.direction() == IteratorDirective.UP)
                || (!before && id.direction() == IteratorDirective.DOWN) ) return reset();
        }
        next = null;
        return null;
    }
    
    /**
     * Ignored.
     */
    public void setAsNeighbor(boolean b) {}
    
    /**
     * Resets iterator to loop from first child to last child of basis.
     */
    public Atom reset() {
        if(basis == null) {next = null; return null;}
        next = basis.node.firstChildAtom();
        last = basis.node.lastChildAtom();
        return next;
    }
    
    /**
     * Returns the next atom in the iteration.
     */
    public Atom next() {
        Atom atom = next;
        next = (next != last) ? next.nextAtom() : null; //check if next is null?
        return atom;
    }        
    
    /**
     * Performs given action for each child atom of basis.
     */
    public void allAtoms(AtomAction act) {
        if(basis == null) return;
        last = basis.node.lastChildAtom();
        for(Atom atom = basis.node.firstChildAtom(); atom != null; atom=atom.nextAtom()) {
            act.actionPerformed(atom);
            if(atom == last) break;
        }
    }
        
    /**
     * Sets the given atom as the basis, so that child atoms of the
     * given atom will be returned upon iteration.  If given atom is
     * a leaf atom, no iterates are given.
     */
    public void setBasis(Atom atom) {
        basis = /*(atom == null || atom.node.isLeaf()) ? null :*/ atom;
    }
    
    /**
     * Returns the current iteration basis.
     */
    public Atom getBasis() {return basis;}
    
    /**
     * The number of atoms returned on a full iteration, using the current basis.
     */
    public int size() {return (basis != null) ? basis.node.childAtomCount() : 0;}   

    private Atom basis;
    private Atom next;
    private Atom last;

}//end of AtomIteratorChildren