package etomica.action;

import etomica.atom.AtomLeaf;
import etomica.atom.IAtom;
import etomica.space.IVector;
import etomica.space.Space;

/**
 * Moves an atom by an amount specified.
 */
public class AtomTranslate extends AtomActionAdapter {
    private static final long serialVersionUID = 1L;
    protected IVector displacement;
        
    public AtomTranslate(Space space) {
        super();
        displacement = space.makeVector();
    }
        
    public final void actionPerformed(IAtom a) {((AtomLeaf)a).getPosition().PE(displacement);}
    public void actionPerformed(IAtom a, IVector d) {((AtomLeaf)a).getPosition().PE(d);}
    public final void setDisplacement(IVector d) {displacement.E(d);}
}