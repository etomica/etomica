package etomica;

/**
 * Iterates over the "neighbors" of a given atom.
 * The neighbors are associated with an atom through entries in its AtomList
 * array.  One entry must be specified for upList neighbors, and another for
 * downList neighbors.
 *
 * @author David Kofke
 */
public class AtomIteratorNeighbor implements AtomIterator {
    
    private IteratorDirective.Direction direction;
    private final IteratorDirective localDirective = new IteratorDirective();
    private final int iUp, iDn;
    private Atom atom;
    private boolean upListNow, doGoDown;
    private final AtomIteratorList iterator;

     /** 
      * Constructor takes indexes to be used to get upNeighbor and downNeighbor
      * lists from an atom.
      */
     public AtomIteratorNeighbor(int iUp, int iDn) {
        iterator = new AtomIteratorList();
        direction = IteratorDirective.UP;
        this.iUp = iUp;
        this.iDn = iDn;
    }
    
    /**
     * Sets the basis for the iterator that loops over the central atoms.
     */
    public void setBasis(Atom a1) {
        atom = a1;
        iterator.reset((AtomLinker)null);//sets hasNext to false
    }
    public Atom getBasis() {return atom;}
    
    /**
     * If setAsNeighbor is true, first returned atom is the one following the reset atom.
     */
    public void setAsNeighbor(boolean b) {iterator.setAsNeighbor(b);}
    
    public final boolean hasNext() {return iterator.hasNext();}
    
    /**
     * 
     */
    public Atom reset(IteratorDirective id) {
        direction = id.direction();
        applyDirection();
        localDirective.copy(id);
        Atom a = null;
        iterator.setBasis((AtomList)null);//set hasNext false
        if(upListNow) {
            iterator.setBasis(atom.atomList[iUp]);
            a = iterator.reset(localDirective.set(IteratorDirective.UP));
        }
        if(!iterator.hasNext() && doGoDown) {
            iterator.setBasis(atom.atomList[iDn]);
            a = iterator.reset(localDirective.set(IteratorDirective.UP));//go up the downlist
        }
        return a;
    }
    
    private void applyDirection() {
        upListNow = direction.doUp();
        doGoDown = direction.doDown();
    }
    
    /**
     * Resets the iterator, so that it is ready to go through all of its iterates,
     * from both the uplist and downlist of the current basis atom.
     */
    public Atom reset() {
        return reset(localDirective.set().set(IteratorDirective.BOTH));
    }
        
    public final Atom next() {
        return nextLinker().atom;
    }
    
    public AtomLinker nextLinker() {
        AtomLinker next = iterator.nextLinker();
        if(!iterator.hasNext() && upListNow && doGoDown) {
            iterator.setBasis(atom.atomList[iDn]);
            iterator.reset(localDirective.set().set(IteratorDirective.UP));//go up the downlist
            upListNow = false;
        }//end if
        return next;
    }
    
    public boolean contains(Atom a) {
        return atom.atomList[iUp].contains(a) || atom.atomList[iDn].contains(a);
    }
    
    public int size() {
        return atom.atomList[iUp].size() + atom.atomList[iDn].size();
    }
    
    /**
     * Not implemented.
     */
    public void allAtoms(AtomAction act) {/* no implementation */}
}  //end of class AtomIteratorNeighbor
    
