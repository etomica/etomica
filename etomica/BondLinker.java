package etomica;

/**
 * Class for making linked-lists of Bond instances.  This is used by the
 * Atoms to list the bonds it is participating in.  Each bond linker
 * points to the bond, the atom using this linker, and linkers preceding
 * and following it in the linked list of bonds for the atom.
 *
 * @author David Kofke
 */
public class BondLinker implements java.io.Serializable {
    
    public Bond bond;
    public BondLinker next, previous;
    public Atom atom;  //that atom that is using this linker
    
    private static BondLinker reservoirFirst = null;
    private static int reservoirCount = 0;
    private static final int reservoirCapacity = 50;
    
    private BondLinker(Bond bond, Atom atom, BondLinker next, BondLinker previous) {
        this.bond = bond;
        this.atom = atom;
        this.next = next;
        this.previous = previous;
        if(next != null) next.previous = this;
        if(previous != null) previous.next = this;
    }
    
    /**
     * Retires this linker (removes it from the list of bonds held
     * by the atom) and repairs the linked list.
     */
    public void delete() {
        if(atom.firstUpBond == this) atom.firstUpBond = next;
        else if(atom.firstDownBond == this) atom.firstDownBond = next;
        
        if(previous != null) previous.next = next;
        if(next != null) next.previous = previous;
        if(reservoirCount < reservoirCapacity) {
            next = reservoirFirst;
            reservoirFirst = this;
            reservoirCount++;
        }
    }
    
    /**
     * Use instead of constructor to obtain a new linker.  Gets one either from
     * a reservoir of unused linkers, or if none available, makes one.
     */
    public static BondLinker getNew(Bond bond, Atom atom, BondLinker next, BondLinker previous) {
        BondLinker newLink = null;
        if(reservoirFirst != null) {//take from reservoir
            newLink = reservoirFirst; //take #1
            reservoirFirst = reservoirFirst.next;//make #2 into #1
            newLink.bond = bond; //fill in fields
            newLink.atom = atom;
            newLink.next = next;
            newLink.previous = previous;
            if(next != null) next.previous = newLink;
            if(previous != null) previous.next = newLink;
            reservoirCount--;
        } else { //reservoir empty; make a new one
            newLink = new BondLinker(bond, atom, next, previous);
        }
        return newLink;
    }
            
}//end of BondLinker
