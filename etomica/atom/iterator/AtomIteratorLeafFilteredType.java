package etomica.atom.iterator;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomList;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IBox;

public class AtomIteratorLeafFilteredType extends AtomIteratorLeafAtoms {

    public AtomIteratorLeafFilteredType(IBox box, IAtomTypeLeaf type) {
        super(box);
        filteredType = type;
    }
    
    public IAtomLeaf nextAtom() {
        IAtomLeaf atom = super.nextAtom();
        while (atom != null) {
            if (atom.getType() == filteredType) {
                return atom;
            }
            atom = super.nextAtom();
        }
        return null;
    }

    public IAtomList next() {
        IAtomList atom = super.next();
        while (atom != null) {
            if (atom.getAtom(0).getType() == filteredType) {
                return atom;
            }
            atom = super.next();
        }
        return null;
    }

    public int size() {
        reset();
        int count = 0;
        for (IAtomLeaf atom = nextAtom(); atom != null; atom = nextAtom()) {
            count++;
        }
        return count;
    }

    private static final long serialVersionUID = 1L;
    protected final IAtomTypeLeaf filteredType;
}
