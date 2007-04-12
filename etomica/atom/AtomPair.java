package etomica.atom;


/**
 * Data structure that contains two mutable atom instances.
 */
public class AtomPair implements AtomSet, java.io.Serializable {

    public AtomPair() {
    }
    
    public AtomPair(IAtom atom0, IAtom atom1) {
        this.atom0 = atom0;
        this.atom1 = atom1;
    }
    
    /* (non-Javadoc)
     * @see etomica.AtomSet#getAtom(int)
     */
    public final IAtom getAtom(int i) {
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

    private static final long serialVersionUID = 1L;
    public IAtom atom0, atom1;

}
