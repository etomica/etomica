package etomica.atom;

import java.io.IOException;

import etomica.space.Space;
import etomica.util.EtomicaObjectInputStream;

 /**
  * Object corresponding to one physical atom or group of atoms. Each atom holds
  * the following publicly accessible fields:
  * <ul>
  * <li>an AtomType instance (fieldname: type) that holds information this atom
  * has in common with other atoms made by the same factory
  * <li>an instance of AtomTreeNode (fieldname: node) that is used to place it
  * in the species hierarchy
  * </ul>
  * <p>
  * @author David Kofke, Andrew Schultz, and C. Daniel Barnes
  * 
  */
public class Atom implements AtomSet, Comparable, java.io.Serializable {

    public Atom(AtomType type, AtomTreeNodeFactory nodeFactory) {
        this.type = type;
        node = nodeFactory.makeNode(this);
    }
    
    /**
     * Makes a simple atom for the given space.  Node is for a leaf atom; 
     * type is a sphere with unit mass and unit size, unique to the new atom; 
     * depth is 0.
     */
    public Atom(Space space) {
        this(new AtomTypeSphere(1,1), AtomTreeNodeLeaf.FACTORY);                        
        node.setIndex(++INSTANCE_COUNT);//default index; changed when added to parent after construction
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
        return type.getAddressManager().sameSpecies(node.getAddress(), atom.node.getAddress());
    }
    /**
     * Returns true if this atom is in the same molecule as the given atom.
     * 
     * @throws NullPointerException
     *             if the argument is null
     */
    public boolean inSameMolecule(Atom atom) {
        return type.getAddressManager().sameMolecule(node.getAddress(), atom.node.getAddress());
    }
    /**
     * Returns true if this atoms is in the same phase as the given atom.
     * 
     * @throws NullPointerException
     *             if the argument is null
     */
    public boolean inSamePhase(Atom atom) {
        return type.getAddressManager().samePhase(node.getAddress(), atom.node.getAddress());
    }
    
    /**
     * Returns a string of digits that uniquely identifies this atom.  String is
     * formed by concatenating the ordinal of this atom to the signature
     * given by the parent of this atom.  If atom has no parent, forms a string
     * from only the ordinal.
     */
    public String signature() {
        if(node.parentGroup() != null) {
            return node.parentGroup().signature() + " " + node.getIndex();
        }
        return Integer.toString(node.getIndex());
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
    	if(node.parentGroup() instanceof SpeciesAgent) return "Molecule(" + signature() + ")";
    	else if(node.isLeaf()) return "Atom(" + signature() + ")";
    	else return "Group(" + signature() + ")";
    }    

    /**
     * Implementation of Comparable interface, with order decided by 
     * placement in Atom hierarchy.  Returns -1, 0, 1 if given Atom
     * is less, equal, or greater, respectively, than this Atom.  
     * Order is determined by compareTo method of atoms' nodes.
     */
    public int compareTo(Object atom) {
        return node.compareTo(((Atom)atom).node);
    }
    
    public void setGlobalIndex(SpeciesMaster speciesMaster) {
        globalIndex = speciesMaster.requestGlobalIndex();
    }

    public final int getGlobalIndex() {
        return globalIndex;
    }
    
    /**
     * Tree node, used to place the atom in the species tree.
     */
    public final AtomTreeNode node;
    
    /**
     * Atom type, holding properties held in common with other atoms made by this atom's
     * factory.
     */
    public final AtomType type;
    
    /**
     * Counter for number of times an atom is instantiated without a parent.  Used
     * to assign a unique index to such atoms.
     */
    private static int INSTANCE_COUNT = 0;
    
    private int globalIndex = -1;
    
    private void readObject(java.io.ObjectInputStream in)
    throws IOException, ClassNotFoundException
    {
        EtomicaObjectInputStream etomicaIn = (EtomicaObjectInputStream)in; 
        etomicaIn.defaultReadObject();
        etomicaIn.addAtom(this);
    }
}
