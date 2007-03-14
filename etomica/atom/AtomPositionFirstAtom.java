package etomica.atom;

import etomica.space.IVector;
import etomica.util.Debug;

/**
 * Returns the position of the first child leaf atom.  Recurses to find
 * the first child leaf atom.
 */

public class AtomPositionFirstAtom implements AtomPositionDefinition, java.io.Serializable {

    public IVector position(Atom atom) {
        AtomLeaf atomLeaf = getFirstChildLeafAtom(atom);
        if (atomLeaf == null) {
            return null;
        }
        return atomLeaf.getCoord().getPosition();
    }
    
    protected AtomLeaf getFirstChildLeafAtom(Atom atom) {
        if (atom.isLeaf()) {
            return (AtomLeaf)atom;
        }
        AtomArrayList childList = ((AtomGroup)atom).getChildList();
        for (int i = 0; i < childList.size(); i++) {
            AtomLeaf a1 = getFirstChildLeafAtom(childList.get(i));
            if(a1 != null) return a1;
            if (Debug.ON) {
                System.out.println("You have yourself a AtomGroup with no children.  That just seems silly.");
            }
        }
        return null;
        
    }

    private static final long serialVersionUID = 1L;
}
