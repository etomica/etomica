package etomica;

/**
* Loops over the atoms in the bond list of a given atom.
*
* @author David Kofke
*/

public class AtomIteratorBonds implements AtomIterator {
    
    private boolean hasNext;
    private boolean upListNow, doGoDown;
    private BondLinker nextBondLink;
    protected IteratorDirective.Direction direction;
    private Atom basis;
    
    public boolean hasNext() {return hasNext;}
    
    public boolean contains(Atom atom) {
        if(basis == null || atom == basis) return false;
        for(BondLinker link=basis.firstUpBond; link!=null; link=link.next) {
            if(link.bond.atom1() == atom || link.bond.atom2() == atom) return true;
        }
        for(BondLinker link=basis.firstDownBond; link!=null; link=link.next) {
            if(link.bond.atom1() == atom || link.bond.atom2() == atom) return true;
        }
        return false;
    }
        
    
    /**
     * Resets iterator to loop up and/or down list of neighbors.  If
     * directive includes an atom specifier, it is used to set the basis
     * for iteration (the basis is the atom whose neighbors are returned by
     * the iterator).
     */
    public Atom reset(IteratorDirective id) {
        direction = id.direction();
        switch(id.atomCount()) {
            case 1: setBasis(id.atom1()); //then fall through to case 0
            case 0: return doReset(); 
            default: hasNext = false;
        }
        return null;
    }
    
    /**
     * Resets iterator to UP using current basis.
     */
    public Atom reset() {
        direction = IteratorDirective.UP;
        return doReset();
    }

    private Atom doReset() {
        upListNow = direction.doUp();
        doGoDown = direction.doDown();
        nextBondLink = null;
        if(basis == null) {
            hasNext = false;
            return null;
        } else {
            if(upListNow) nextBondLink = basis.firstUpBond;
            if(nextBondLink == null && doGoDown) {
                nextBondLink = basis.firstDownBond;
                upListNow = false;
            }
        }
        hasNext = (nextBondLink != null);
        if(hasNext) return (nextBondLink.bond.link1==nextBondLink) ? nextBondLink.bond.link2.atom : nextBondLink.bond.link1.atom;
        else return null;
    }//end doReset
        
    /**
     * Has no effect, since this iterator always behaves as a neighbor iterator.
     */
    public void setAsNeighbor(boolean b) {
        //need an exception for this
//        if(!b) System.out.println("Error in AtomIteratorBonds; attempt to set a not a neighborIterator");
    }
    
    public Atom next() {
        BondLinker next = nextBondLink;
        nextBondLink = nextBondLink.next;
        hasNext = nextBondLink != null;
        if(!hasNext && upListNow && doGoDown) {//done going up and now prepare to go down
            nextBondLink = basis.firstDownBond;
            hasNext = (nextBondLink != null);
            upListNow = false;
        }
        return (next.bond.link1==next) ? next.bond.link2.atom : next.bond.link1.atom;//need instead to loop over atoms of bond
    }
    
    //not implemented
    public void allAtoms(AtomAction act) {
        System.out.println("allAtoms not implemented in AtomIteratorBonds");
        System.exit(1);
    }

    /**
     * Sets basis to given atom, and sets hasNext to false (must perform a separate call to reset
     * before commencing iteration).
     */
    public void setBasis(Atom atom) {
        basis = atom;
        hasNext = false;
    }

}//end of AtomIteratorBonds