package etomica.action;

import etomica.api.IAtom;
import etomica.api.IAtomPositioned;
import etomica.api.IVector;
import etomica.space.Space;

/**
 * Moves an atom by an amount specified.
 */
public class AtomTranslate implements AtomAction {
    private static final long serialVersionUID = 1L;
    protected IVector displacement;
        
    public AtomTranslate(Space space) {
        super();
        displacement = space.makeVector();
    }
        
    public final void actionPerformed(IAtom a) {((IAtomPositioned)a).getPosition().PE(displacement);}
    public void actionPerformed(IAtom a, IVector d) {((IAtomPositioned)a).getPosition().PE(d);}
    public final void setDisplacement(IVector d) {displacement.E(d);}
}