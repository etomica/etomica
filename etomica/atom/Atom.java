package etomica.atom;

import java.io.IOException;

import etomica.chem.elements.ElementSimple;
import etomica.phase.Phase;
import etomica.util.Debug;
import etomica.util.EtomicaObjectInputStream;

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
public abstract class Atom implements AtomSet, Comparable, java.io.Serializable {

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
     * Returns 1, indicating that this AtomSet is an Atom.
     */
    public final int count() {return 1;}
    
    /**
     * Returns this if i==0, otherwise throws exception.
     * 
     * @throws IllegalArgumentException if i != 0
     */
    public final Atom getAtom(int i) {
        if(i == 0 ) return this;
        throw new IllegalArgumentException();
    }
    
    /**
     * Returns true if the given object is this atom instance, or if it is
     * a length-1 AtomSet holding this instance.
     */
    public boolean equals(Object object) {
        if(!(object instanceof AtomSet) || ((AtomSet)object).count() != 1) return false;
        return this == ((AtomSet)object).getAtom(0);
    }

    /**
     * Returns true if the given object is this atom instance, or if it is
     * a length-1 AtomSet holding this instance.
     */
    public final boolean equals(AtomSet atoms) {
        if (atoms == null || atoms.count() != 1) return false;
        return this == atoms.getAtom(0);
    }

    /**
     * Returns true if this atom is in the same species as the given atom.
     * 
     * @throws NullPointerException
     *             if the argument is null
     */
    public boolean inSameSpecies(Atom atom) {
        return type.getAddressManager().sameSpecies(getAddress(), atom.getAddress());
    }
    
    /**
     * Returns true if this atom is in the same molecule as the given atom.
     * 
     * @throws NullPointerException
     *             if the argument is null
     */
    public boolean inSameMolecule(Atom atom) {
        return type.getAddressManager().sameMolecule(getAddress(), atom.getAddress());
    }
    /**
     * Returns true if this atoms is in the same phase as the given atom.
     * 
     * @throws NullPointerException
     *             if the argument is null
     */
    public boolean inSamePhase(Atom atom) {
        return type.getAddressManager().samePhase(getAddress(), atom.getAddress());
    }
    
    /**
     * Returns a string of digits that uniquely identifies this atom.  String is
     * formed by concatenating the ordinal of this atom to the signature
     * given by the parent of this atom.  If atom has no parent, forms a string
     * from only the ordinal.
     */
    public String signature() {
        if(parent != null) {
            return parent.signature() + " " + getIndex();
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
    	if(this instanceof SpeciesMaster) return "Master(" + signature() + ")";
    	else if(this instanceof SpeciesAgent) return "Agent(" + signature() + ")";
    	if(parent instanceof SpeciesAgent) return "Molecule(" + signature() + ")";
    	else if(isLeaf()) return "Atom(" + signature() + ")";
    	else return "Group(" + signature() + ")";
    }    

    /**
     * Implementation of Comparable interface, with order decided by 
     * placement in Atom hierarchy.  Returns -1, 0, 1 if given Atom
     * is less, equal, or greater, respectively, than this Atom.  
     * Order is determined by compareTo method of atoms' nodes.
     */
    public int compareTo(Object otherAtom) {
        int otherAddress = ((Atom)otherAtom).atomTreeAddress;
        //use "<" for ">" because atomIndex is negative
        return otherAddress < atomTreeAddress ? 1 : (otherAddress == atomTreeAddress ? 0 : -1);
    }
    
    public void setGlobalIndex(SpeciesMaster speciesMaster) {
        globalIndex = speciesMaster.requestGlobalIndex();
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

    public abstract boolean isLeaf();
    
    public void dispose() {
        setParent(null);
    }

    public void setParent(AtomGroup newParent) {
        
        //old parent is not null; remove from it
        if(parent != null) {
            int index = getIndex();
            parent.getChildList().removeAndReplace(index);
            parent.getChildList().maybeTrimToSize();
            if (parent.getChildList().size() > index) {
                parent.getChildList().get(index).setIndex(index);
            }
            parent.removeAtomNotify(this);
        }
        
        parent = newParent;

        if(parent == null) {//new parent is null
            return;
        }

        setIndex(parent.getChildList().size());
        
        parent.getChildList().add(this);
        parent.addAtomNotify(this);
    }//end of addAtom

    public AtomGroup getParentGroup() {
        return parent;
    }
    
    /**
     * Returns the molecule in which this atom resides.  A "molecule" is an atomgroup
     * that is one step below a species agent in the hierarchy of atomgroups.
     */
    public Atom getParentMolecule() {
        if(parent == null) return null;
        return (parent instanceof SpeciesAgent) ? this : parent.getParentMolecule();
    }
                
    /**
     * Phase in which this atom resides
     */
    public Phase getParentPhase() {
        return (parent != null) ? parent.getParentPhase() : null;
    }
    
    public SpeciesAgent getParentSpeciesAgent() {
        return (parent != null) ? parent.getParentSpeciesAgent() : null;
    }
    
    /**
     * Integer assigned to this atom by its parent molecule.
     */
    public final int getAddress() {return atomTreeAddress;}

    public void setIndex(int index) {
        setIndex((parent != null) ? parent.getAddress() : 0, index);
    }
    
    protected void setIndex(int parentAddress, int index) {
        atomTreeAddress = parentAddress + type.getAddressManager().shiftIndex(index);
        if (Debug.ON && type.getAddressManager().getIndex(atomTreeAddress) != index) {
            type.getAddressManager().getIndex(atomTreeAddress);
            throw new RuntimeException(atomTreeAddress+" "+index+" "+(type.getAddressManager().getIndex(atomTreeAddress)));
        }
    }
    
    public int getIndex() {
        return type.getAddressManager().getIndex(atomTreeAddress);
    }
    
    public int getMoleculeIndex() {
        return type.getAddressManager().getMoleculeIndex(atomTreeAddress);
    }

    public int getPhaseIndex() {
        return type.getAddressManager().getPhaseIndex(atomTreeAddress);
    }

    public int getSpeciesIndex() {
        return type.getAddressManager().getSpeciesIndex(atomTreeAddress);
    }

    public boolean inSimulation() {
        return atomTreeAddress < 0;
    }
    
    /**
     * Returns true if the given atom is in the hierarchy of parents of this atom,
     * or if the given atom is this atom.  Returns true, for example, if the given
     * atom is this atom's parent, or its parent's parent, etc.
     */ 
    public boolean isDescendedFrom(Atom group) {
        if(group == null) return false;
        return group.getType().getAddressManager().sameAncestry(group.atomTreeAddress,this.atomTreeAddress);
    }
    
    /**
     * Returns the child of the given node through which this node is derived.
     * If given node is parent node of this, returns this.
     * If this node is not descended from the given node, returns null.
     */
    public Atom getChildWhereDescendedFrom(Atom atom) {
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
    protected AtomGroup parent;
    
    private void readObject(java.io.ObjectInputStream in)
    throws IOException, ClassNotFoundException
    {
        EtomicaObjectInputStream etomicaIn = (EtomicaObjectInputStream)in; 
        etomicaIn.defaultReadObject();
        etomicaIn.addAtom(this);
    }
}
