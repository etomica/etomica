package etomica.atom;

import etomica.space.IVector;
import etomica.space.Space;
import etomica.util.Debug;

 /**
  * Atom that represents a physical atom with a position.
  * <p>
  * @author David Kofke, Andrew Schultz, and C. Daniel Barnes
  * 
  */
public class AtomLeaf extends Atom implements IAtomLeaf, IAtomPositioned {

    public AtomLeaf(Space space, AtomTypeLeaf type) {
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

//    public void setIndex(int index) {
//        setIndex((parent != null) ? parent.getAddress() : 0, index);
//    }
//    
    /**
     * Returns a string of digits that uniquely identifies this atom.  String is
     * formed by concatenating the ordinal of this atom to the signature
     * given by the parent of this atom.  If atom has no parent, forms a string
     * from only the ordinal.
     */
    public String signature() {
        if(parent != null) {
            if (parent instanceof Atom) {
                return ((Atom)parent).signature() + " " + getIndex();
            }
            return parent.toString() + " " + getIndex();
        }
        return Integer.toString(getIndex());
    }
    
    /**
     * Returns a string formed by concatenating the signature of this atom
     * to a string that identifies it as a species master, species agent, 
     * molecule, group, or (leaf) atom.
     */
    public final String toString() {
//        return Integer.toBinaryString(node.index());
        return "Atom(" + signature() + ")";
    }

    /**
     * Informs the Atom that the given AtomGroup is its parent.
     * This method should only be called by the parent.
     */
    public void setParent(IMolecule newParent) {
        if (Debug.ON && ((newParent == null) == (parent == null))) {
            // newParent must be null (removal) if current parent is not null
            // new Parent must not be null if current parent null
            throw new IllegalArgumentException(newParent+" is not my parent");
        }
        parent = newParent;
    }

    public IMolecule getParentGroup() {
        return parent;
    }
    
    public IVector getPosition() {
        return position;
    }
    
    private static final long serialVersionUID = 2L;
    protected final IVector position;
    protected IMolecule parent;
}
