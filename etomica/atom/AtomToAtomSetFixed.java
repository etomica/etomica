package etomica.atom;

import java.io.Serializable;

import etomica.api.IAtom;
import etomica.api.IAtomList;

public class AtomToAtomSetFixed implements AtomToAtomSet, AtomToIndex, Serializable {

    private static final long serialVersionUID = 1L;

    public AtomToAtomSetFixed() {
        atomArrayList = new AtomArrayList();
    }
    
    public void setArrayList(AtomArrayList list) {
        if (list == null) {
            atomArrayList = new AtomArrayList();
        }
        else {
            atomArrayList = list;
        }
    }
    
    public IAtomList getAtomSet(IAtom atom) {
        return atomArrayList;
    }
    
    public int getIndex(IAtom atom) {
        return atomArrayList.indexOf(atom);
    }

    private AtomArrayList atomArrayList;
}
