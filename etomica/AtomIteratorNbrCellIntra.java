package etomica;

/**
* Iterates among the children of a given basis, those atoms
* that are cell-list neighbors of a specified atom that is
* a child of the same basis.
*
* @author David Kofke
*/
//would like to modify so that central atom can be any descendant of the basis.
public final class AtomIteratorNbrCellIntra implements AtomIterator {
    
    /**
     * Indicates if another iterate is forthcoming.
     */
    public boolean hasNext() {return nextAtom != null;}
    
    /**
     * True if the parent group of the given atom is the current basis for the iterator.
     * False otherwise, or if atom or basis is null.
     */
    public boolean contains(Atom atom) {
        return atom != null && basis != null && atom.node.parentGroup() == basis;
    }
    
    /**
     * Does reset if atom in iterator directive is child of the current basis.  
     * Sets hasNext false if given atom does is not child of basis.  Throws
     * an IllegalArgumentException if directive does not specify an atom.
     */
    public Atom reset(IteratorDirective id) {
        direction = id.direction();
        return reset(id.atom1());
    }
    
    public Atom reset(Atom atom) {
        referenceAtom = atom;
        upListNow = direction.doUp();
        doGoDown = direction.doDown();
        nextAtom = null;
        if(atom == null) {
            throw new IllegalArgumentException("Cannot reset AtomIteratorNeighborIntra without referencing an atom");
        //probably need isDescendedFrom instead of parentGroup here
        } else if(atom.node.parentGroup() == basis) {
            cell = ((AtomSequencerCell)atom.seq).site();
//                cellIterator = nextCoordinate.cell.neighborIterator(); //cell-neighbor iterator for cell containing this atom
//                cellIterator.resetUp();  //next cell returned by iterator is first up the neighbor list
            if(upListNow) {
                cellIterator.reset(cell, IteratorDirective.UP);//set cell iterator to return next cell up (shouldn't begin with this cell)
                nextAtom = atom.seq.nextAtom();
            }
            if(nextAtom == null) advanceCell();
        }
        return nextAtom;
    }
                
    // Finds first atom of next occupied cell
    private void advanceCell() {
        do {
            if(cellIterator.hasNext()) {
                nextAtom = upListNow ? ((AtomCell)cellIterator.next()).first()
                                     : ((AtomCell)cellIterator.next()).last();
            } else if(doGoDown) {//no more cells that way; see if should now reset to look at down-cells
                cellIterator.reset(cell, IteratorDirective.DOWN);//set cell iterator to return next cell down
                nextAtom = referenceAtom.seq.previousAtom();
                upListNow = false;
                doGoDown = false;
            } else {//no more cells at all
                break;
            }
        } while(nextAtom == null);
    }
            
    public Atom next() {
        Atom atom = nextAtom;
        nextAtom = upListNow ? atom.seq.nextAtom() : atom.seq.previousAtom();
        if(nextAtom == null) advanceCell();
        return atom;
    }
    /**
     * Ignored.
     */
    public void setAsNeighbor(boolean b) {}
    
    /**
     * Resets iterator to loop from first child to last child of basis.
     */
    public Atom reset() {
        throw new IllegalArgumentException("Cannot reset AtomIteratorNeighborIntra without referencing an atom");
        if(basis == null) {next = null; return null;}
        next = basis.node.firstChildAtom();
        last = basis.node.lastChildAtom();
        return next;
    }
    
    
    /**
     * Performs given action for each child atom of basis.
     */
    public void allAtoms(AtomAction act) {
        throw new RuntimeException("AtomIteratorNbrCellIntra.allAtoms not implemented");
/*        if(basis == null) return;
        last = basis.node.lastChildAtom();
        for(Atom atom = basis.node.firstChildAtom(); atom != null; atom=atom.nextAtom()) {
            act.actionPerformed(atom);
            if(atom == last) break;
        }*/
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
    private Atom referenceAtom;
    private boolean upListNow, doGoDown;
    private IteratorDirective.Direction direction;
    private AtomCell cell;

}//end of AtomIteratorNbrCellIntra