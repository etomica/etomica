package etomica.atom;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomList;


/**
 * Data structure that contains two mutable atom instances.
 */
public class AtomPair implements IAtomList, java.io.Serializable {

    public AtomPair() {
    }
    
    public AtomPair(IAtomLeaf atom0, IAtomLeaf atom1) {
        this.atom0 = atom0;
        this.atom1 = atom1;
    }
    
    /* (non-Javadoc)
     * @see etomica.AtomSet#getAtom(int)
     */
    public final IAtomLeaf getAtom(int i) {
        if(i == 0) return atom0;
        if(i == 1) return atom1;
        throw new IllegalArgumentException();
    }

    /* (non-Javadoc)
     * @see etomica.AtomSet#count()
     */
    public final int getAtomCount() {
        return 2;
    }
    
    public String toString() {
        return "["+atom0+","+atom1+"]";
    }

    private static final long serialVersionUID = 1L;
    public IAtomLeaf atom0, atom1;

}
