package etomica.atom;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IMolecule;
import etomica.api.IVectorMutable;
import etomica.space.ISpace;
import etomica.util.Debug;

 /**
  * Atom that represents a physical atom with a position.
  * <p>
  * @author David Kofke, Andrew Schultz, and C. Daniel Barnes
  * 
  */
public class AtomLeaf extends Atom implements IAtomLeaf, IAtomPositioned {

    public AtomLeaf(ISpace space, IAtomTypeLeaf type) {
        super();
        this.type = type;
        position = space.makeVector();
    }
    
    /**
     * Makes a simple atom for the given space.  Coordinate is non-kinetic sphere;
     * node is for a leaf atom; type is a sphere with unit mass and unit size, 
     * unique to the new atom; depth is 0.
     */
    public AtomLeaf(ISpace space) {
        super();
        type = null;
        position = space.makeVector();
    }

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
        return "Atom(" + signature() + ")";
    }

    public final void setIndex(int newIndex) {
        index = newIndex;
    }
    
    public final int getIndex() {
        return index;
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
    
    public IVectorMutable getPosition() {
        return position;
    }
    
    public int getLeafIndex() {
        return leafIndex;
    }
    
    public void setLeafIndex(int newLeafIndex) {
        leafIndex = newLeafIndex;
    }

    /**
     * @return the Atom type, holding properties held in common with other 
     * atoms made by this atom's factory.
     */
    public final IAtomTypeLeaf getType() {
        return type;
    }

    private static final long serialVersionUID = 3L;
    protected final IAtomTypeLeaf type;
    protected int index;
    protected final IVectorMutable position;
    protected IMolecule parent;
    protected int leafIndex;
}
