package etomica.atom;


/**
 * Data structure that contains two mutable atom instances.
 */
public class AtomPair implements AtomSet, Comparable, java.io.Serializable {

    public AtomPair() {
    }
    
    public AtomPair(Atom atom0, Atom atom1) {
        this.atom0 = atom0;
        this.atom1 = atom1;
    }
    
    /* (non-Javadoc)
     * @see etomica.AtomSet#getAtom(int)
     */
    public final Atom getAtom(int i) {
        if(i == 0) return atom0;
        if(i == 1) return atom1;
        throw new IllegalArgumentException();
    }

    /* (non-Javadoc)
     * @see etomica.AtomSet#count()
     */
    public final int count() {
        return 2;
    }
    
    public void copyTo(AtomPair pair) {
        pair.atom0 = atom0;
        pair.atom1 = atom1;
    }
    
    public String toString() {
        return "["+atom0.toString()+","+atom1.toString()+"]";
    }

    /**
     * Returns result of compareTo applied between first atoms 
     * of the two pairs, and if that is zero, returns comparison
     * of second atoms of pair.  Implementation of Comparable interface.
     */
    public int compareTo(Object pair) {
        int i0 = atom0.compareTo(((AtomPair)pair).atom0);
        return (i0 != 0) ? i0 : atom1.compareTo(((AtomPair)pair).atom1);
    }
    
    private static final long serialVersionUID = 1L;
    public Atom atom0, atom1;

}
