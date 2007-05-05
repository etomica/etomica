package etomica.atom;

import etomica.space.IVector;
import etomica.space.Space;

 /**
  * Atom that represents a physical atom with a position.
  * <p>
  * @author David Kofke, Andrew Schultz, and C. Daniel Barnes
  * 
  */
public class AtomLeaf extends Atom implements IAtomPositioned {

    public AtomLeaf(Space space, AtomType type) {
        super(type);
        position = space.makeVector();
    }
    
    /**
     * Makes a simple atom for the given space.  Coordinate is non-kinetic sphere;
     * node is for a leaf atom; type is a sphere with unit mass and unit size, 
     * unique to the new atom; depth is 0.
     */
    public AtomLeaf(Space space) {
        super();
        position = space.makeVector();
    }

    public IVector getPosition() {
        return position;
    }
    
    private static final long serialVersionUID = 2L;
    protected final IVector position;
}
