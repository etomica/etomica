package etomica.action;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Moves an atom by an amount specified.
 */
public class AtomTranslate extends AtomActionAdapter {
    private static final long serialVersionUID = 1L;
    protected Vector displacement;
        
    public AtomTranslate(Space space) {
        super();
        displacement = space.makeVector();
    }
        
    public final void actionPerformed(Atom a) {((AtomLeaf)a).getCoord().position().PE(displacement);}
    public void actionPerformed(Atom a, Vector d) {((AtomLeaf)a).getCoord().position().PE(d);}
    public final void setDisplacement(Vector d) {displacement.E(d);}
}