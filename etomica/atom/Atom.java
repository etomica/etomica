package etomica.atom;

import etomica.chem.elements.ElementSimple;
import etomica.util.Debug;

 /**
  * Object corresponding to one physical atom or group of atoms. Each atom holds
  * the following publicly accessible fields:
  * <ul>
  * <li>an AtomType instance (fieldname: type) that holds information this atom
  * has in common with other atoms made by the same factory
  * </ul>
  * <p>
  * @author David Kofke, Andrew Schultz, and C. Daniel Barnes
  * 
  */
public abstract class Atom implements IAtom, java.io.Serializable {

    public Atom(AtomType type) {
        this.type = type;
    }
    
    /**
     * Makes a simple atom.  Node is for a leaf atom; 
     * type is a sphere with unit mass and unit size, unique to the new atom; 
     * depth is 0.
     */
    public Atom() {
        this(makeAtomTypeSphere());                        
        setIndex(++INSTANCE_COUNT);//default index; changed when added to parent after construction
    }
    
    /**
     * Method to return a dummy AtomType that's valid (has parent explicitly set to null)
     */
    private static AtomTypeSphere makeAtomTypeSphere() {
        AtomTypeSphere newType = new AtomTypeSphere(new ElementSimple("Simple",1),1);
        newType.setParentType(null);
        return newType;
    }

    /**
     * Returns true if this atom is in the same molecule as the given atom.
     * 
     * @throws NullPointerException
     *             if the argument is null
     */
    public boolean inSameMolecule(IAtom atom) {
        return type.getAddressManager().sameMolecule(getAddress(), atom.getAddress());
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
//        return Integer.toBinaryString(node.index());
    	if(this instanceof ISpeciesAgent) return "Agent(" + signature() + ")";
    	if(parent instanceof ISpeciesAgent) return "Molecule(" + signature() + ")";
    	else if(!(this instanceof IAtomGroup)) return "Atom(" + signature() + ")";
    	else return "Group(" + signature() + ")";
    }    

    public void setGlobalIndex(AtomManager atomManager) {
        globalIndex = atomManager.requestGlobalIndex();
    }

    public final int getGlobalIndex() {
        return globalIndex;
    }
    
    /**
     * @return the Atom type, holding properties held in common with other 
     * atoms made by this atom's factory.
     */
    public final AtomType getType() {
        return type;
    }

    /**
     * Informs the Atom that the given AtomGroup is its parent.
     * This method should only be called by the parent.
     */
    public void setParent(IAtomGroup newParent) {
        if (Debug.ON && ((newParent == null) == (parent == null))) {
            // newParent must be null (removal) if current parent is not null
            // new Parent must not be null if current parent null
            throw new IllegalArgumentException(newParent+" is not my parent");
        }
        parent = newParent;
    }

    public IAtomGroup getParentGroup() {
        return parent;
    }
    
    /**
     * Integer assigned to this atom by its parent molecule.
     */
    public final int getAddress() {return atomTreeAddress;}

    public void setIndex(int index) {
        setIndex((parent != null) ? parent.getAddress() : 0, index);
    }
    
    private void setIndex(int parentAddress, int index) {
        atomTreeAddress = parentAddress + type.getAddressManager().shiftIndex(index);
        if (Debug.ON && type.getAddressManager().getIndex(atomTreeAddress) != index) {
            type.getAddressManager().getIndex(atomTreeAddress);
            throw new RuntimeException(atomTreeAddress+" "+index+" "+(type.getAddressManager().getIndex(atomTreeAddress)));
        }
    }
    
    public int getIndex() {
        return type.getAddressManager().getIndex(atomTreeAddress);
    }
    
    /**
     * Returns true if the given atom is in the hierarchy of parents of this atom,
     * or if the given atom is this atom.  Returns true, for example, if the given
     * atom is this atom's parent, or its parent's parent, etc.
     */ 
    public boolean isDescendedFrom(IAtom group) {
        if(group == null) return false;
        return group.getType().getAddressManager().sameAncestry(group.getAddress(),atomTreeAddress);
    }
    
    /**
     * Returns the child of the given node through which this node is derived.
     * If given node is parent node of this, returns this.
     * If this node is not descended from the given node, returns null.
     */
    public IAtom getChildWhereDescendedFrom(IAtom atom) {
        if(parent == null) return null;
        return (parent == atom) ? this : parent.getChildWhereDescendedFrom(atom);
    }
    
    protected final AtomType type;
    
    /**
     * Counter for number of times an atom is instantiated without a parent.  Used
     * to assign a unique index to such atoms.
     */
    private static int INSTANCE_COUNT = 0;
    
    private int globalIndex = -1;
    protected int atomTreeAddress;
    protected IAtomGroup parent;
}
