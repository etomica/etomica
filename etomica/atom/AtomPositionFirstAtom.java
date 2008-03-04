package etomica.atom;

import etomica.api.IAtom;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IMolecule;
import etomica.api.IVector;
import etomica.util.Debug;

/**
 * Returns the position of the first child leaf atom.  Recurses to find
 * the first child leaf atom.
 */

public class AtomPositionFirstAtom implements AtomPositionDefinition, java.io.Serializable {

    public IVector position(IAtom atom) {
        IAtomPositioned atomLeaf = getFirstChildLeafAtom(atom);
        if (atomLeaf == null) {
            return null;
        }
        return atomLeaf.getPosition();
    }
    
    protected IAtomPositioned getFirstChildLeafAtom(IAtom atom) {
        if (!(atom instanceof IMolecule)) {
            return (IAtomPositioned)atom;
        }
        IAtomSet childList = ((IMolecule)atom).getChildList();
        for (int i = 0; i < childList.getAtomCount(); i++) {
            IAtomPositioned a1 = getFirstChildLeafAtom(childList.getAtom(i));
            if(a1 != null) return a1;
            if (Debug.ON) {
                System.out.println("You have yourself a AtomGroup with no children.  That just seems silly.");
            }
        }
        return null;
        
    }

    private static final long serialVersionUID = 1L;
}
