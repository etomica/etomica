package etomica.atom;

import java.io.Serializable;

public class AtomToArrayListFixed implements AtomToArrayList, AtomToIndex, Serializable {

    private static final long serialVersionUID = 1L;

    public AtomToArrayListFixed() {
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
    
    public AtomArrayList getArrayList(IAtom atom) {
        return atomArrayList;
    }
    
    public int getIndex(IAtom atom) {
        return atomArrayList.indexOf(atom);
    }

    private AtomArrayList atomArrayList;
}
