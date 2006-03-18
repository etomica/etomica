package etomica.atom;

import java.io.Serializable;

public class AtomToArrayListFixed implements AtomToArrayList, AtomToIndex, Serializable {

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
    
    public AtomArrayList getArrayList(Atom atom) {
        return atomArrayList;
    }
    
    public int getIndex(Atom atom) {
        return atomArrayList.indexOf(atom);
    }

    private AtomArrayList atomArrayList;
}
