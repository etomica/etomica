package etomica.atom;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomList;


/**
 * Data structure that contains a single mutable atom instance.
 */
public class AtomSetSinglet implements IAtomList, java.io.Serializable {

    public AtomSetSinglet() {
    }
    
    public AtomSetSinglet(IAtomLeaf atom) {
        this.atom = atom;
    }
    
    public final IAtomLeaf getAtom(int i) {
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
    public IAtomLeaf atom;
}
