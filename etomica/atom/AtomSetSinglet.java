package etomica.atom;


/**
 * Data structure that contains a single mutable atom instance.
 */
public class AtomSetSinglet implements AtomSet, java.io.Serializable {

    public AtomSetSinglet() {
    }
    
    public AtomSetSinglet(IAtom atom) {
        this.atom = atom;
    }
    
    public final IAtom getAtom(int i) {
        if(i == 0) return atom;
        throw new IllegalArgumentException();
    }

    public final int getAtomCount() {
        return 1;
    }
    
    public String toString() {
        return "["+atom+"]";
    }

    private static final long serialVersionUID = 1L;
    public IAtom atom;
}
