package etomica.atom;

import etomica.species.Species;
import etomica.util.Debug;
import etomica.util.Default;

//import etomica.electrostatics.*;

/**
 * AtomType holds fields that are common to many atoms. It serves many
 * functions:
 * <ul>
 * <li>it holds parameters that define the general features of the atom (e.g.,
 * size).
 * <li>it holds an index manager that is used to interpret the atom's index
 * held in its node
 * </ul>
 * The AtomType of an atom is set by its AtomFactory when it builds a molecule.
 * Each Atom has an instance variable named "type" that holds the AtomType
 * instance.
 * <p>
 * AtomType instances are arranged in a hierarchy that parallels the species
 * hierarchy. The AtomType of an atom group will be an instance of
 * AtomTypeGroup, and will have a set of child types that are the types of the
 * child atoms of the atom group. Whereas an atom group may have more than one
 * child of a given type, an AtomTypeGroup will have the AtomType for those
 * atoms represented only once.
 */

public abstract class AtomType implements java.io.Serializable, Comparable {

    AtomFactory creator;//set in constructor of AtomFactory
    protected int speciesIndex = -1;
    private Species species;
    private int index;
    private int typeTreeAddress;

    private AtomAddressManager addressManager;
    private AtomPositionDefinition positionDefinition;

    protected AtomTypeGroup parentType;
    
    private boolean isInteracting = false;

    /**
     * Used only to create root type.
     */
    AtomType(AtomAddressManager indexManager) {
        parentType = null;
        this.addressManager = indexManager;
        setChildIndex(0,0);
        positionDefinition = null;
        index = 0;
    }

    /**
     * Primary constructor called by subclasses to make the instance.
     * 
     * @param parentType
     *            the instance of the AtomType that is the parent of this in the
     *            AtomType hierarchy. This constructor adds this instance to the
     *            list of child types of the given parent (if not null, which is
     *            permitted and may arise if working outside the species tree
     *            hierarchy)
     * @param positionDefinition
     *            used by many classes as default choice for defining the
     *            spatial position of an instance of an atom(group) of this type
     */
    public AtomType(AtomPositionDefinition positionDefinition) {
        this.positionDefinition = positionDefinition;
        setParentType(null);
    }

    /**
     * Integer describing this AtomType's place in the tree
     */
    public final int getAddress() {
        return typeTreeAddress;
    }

    /**
     * @param parentType
     *            the instance of the AtomType that is the parent of this in the
     *            AtomType hierarchy. This constructor adds this instance to the
     *            list of child types of the given parent (if not null, which is
     *            permitted and may arise if working outside the species tree
     *            hierarchy)
     */
    public void setParentType(AtomTypeGroup newParentType) {
        if (parentType != null) {
            throw new RuntimeException("You can't set the parent type twice");
        }
        
        parentType = newParentType;
        if (parentType == null) {
            addressManager = AtomAddressManager
                    .makeSimpleIndexManager(Default.BIT_LENGTH);
            index = 0;
            setChildIndex(0,0);
        } else {
            addressManager = parentType.getAddressManager().makeChildManager();
            setChildIndex(parentType.childTypes.length);
            index = parentType.requestIndex();
            parentType.addChildType(this);
        }
    }        

    public void setChildIndex(int newChildIndex) {
        setChildIndex((parentType != null) ? parentType.getAddress() : 0, newChildIndex);
    }
    
    protected void setChildIndex(int parentTypeAddress, int newChildIndex) {
        typeTreeAddress = parentTypeAddress + addressManager.shiftIndex(newChildIndex);
        if (Debug.ON && getAddressManager().getIndex(typeTreeAddress) != newChildIndex) {
            addressManager.getIndex(typeTreeAddress);
            throw new RuntimeException(typeTreeAddress+" "+newChildIndex+" "+(addressManager.getIndex(typeTreeAddress)));
        }
    }
    
    public int getChildIndex() {
        return addressManager.getIndex(typeTreeAddress);
    }
    
    public int getIndex() {
        return index;
    }
    
    void resetIndex() {
        index = parentType.requestIndex();
    }
    
    /**
     * Indicates whether this type is at the bottom of the AtomType hierarchy.
     * 
     * @return true if this type is an instance of AtomTypeLeaf
     */
    public abstract boolean isLeaf();

    /**
     * Returns the parent AtomType of this AtomType.
     */
    public AtomTypeGroup getParentType() {
        return parentType;
    }

    /**
     * Returns the index manager used to set and interpret the 
     * index of an atom of this type.
     */
    public AtomAddressManager getAddressManager() {
        return addressManager;
    }

    /**
     * Implementation of Comparable interface, based on type index value.  
     * Returns -1, 0, 1 if getIndexManager().getTypeIndex() for this type
     * is less, equal, or greater, respectively, than the corresponding 
     * value for the given type.
     */
    public int compareTo(Object atomType) {
        int otherAddress = ((AtomType) atomType).getAddress();
        return otherAddress > typeTreeAddress ? -1 : (otherAddress == typeTreeAddress ? 0 : 1);
    }

    /**
     * The position definition held by the type provides an appropriate default
     * to define the position of an atom of this type. This field is set in the
     * definition of the parent species of the atom. It is null for SpeciesRoot,
     * SpeciesMaster, and SpeciesAgent atoms.
     * 
     * @return Returns the PositionDefinition for an atom of this type.
     */
    public AtomPositionDefinition getPositionDefinition() {
        return positionDefinition;
    }

    /**
     * Sets the PositionDefinition used for this AtomType
     */
    public void setPositionDefinition(AtomPositionDefinition newPositionDefinition) {
        positionDefinition = newPositionDefinition;
    }
    
    /**
     * Returns true if an atom of this type is descended from an atom (any atom)
     * having the given type.
     */
    public boolean isDescendedFrom(AtomType type) {
        return type.getAddressManager().sameAncestry(type.getAddress(),typeTreeAddress);
    }

    /**
     * Returns the AtomFactory that creates an atom of this type.
     */
    public AtomFactory creator() {
        return creator;
    }

    /**
     * Returns the depth of this atom in the atom hierarchy. That is, returns
     * the number of parent relations between this atom and the species master.
     */
    public final int getDepth() {
        return addressManager.getDepth();
    }

    /**
     * @return Returns the species.
     */
    public Species getSpecies() {
        return species;
    }

    /**
     * @param species
     *            The species to set.
     */
    public void setSpecies(Species species) {
        this.species = species;
        speciesIndex = species.getIndex();
    }

    /**
     * @return Returns the speciesIndex.
     */
    public final int getSpeciesIndex() {
        return speciesIndex;
    }

    public void setInteracting(boolean b) {
        isInteracting = b;
    }

    /**
     * Returns true if one or more potentials are defined to act
     * on an atom of this type.
     */
    public boolean isInteracting() {
        return isInteracting;
    }
    
    public String toString() {
        return "AtomType "+index;
    }
    
    //prototype of a real atom type
    /*
     * public final static class Carbon extends Sphere { public Carbon() {
     * super(12.0, Color.black, 1.1); //mass, color, diameter } } public final
     * static class Carbon12 extends Sphere { public Carbon12() { super(12.0,
     * Color.black, 1.1); this.setName("Carbon" + Integer.toString(CarbonID++)); } }
     * 
     * public final static class Hydrogen extends Sphere { public Hydrogen() {
     * super(1.0, Color.cyan, 0.5); this.setName("Hydrogen" +
     * Integer.toString(HydrogenID++)); } }
     * 
     * public final static class Oxygen extends Sphere { public Oxygen() {
     * super(16.0, Color.red, 1.3); this.setName("Oxygen" +
     * Integer.toString(OxygenID++)); } }
     * 
     * public final static class Nitrogen extends Sphere { public Nitrogen() {
     * super(14.0, Color.blue, 1.2); this.setName("Nitrogen" +
     * Integer.toString(NitrogenID++)); } }
     *  
     */
    //interfaces for anisotropic atom types
    public interface Rotator {

        public double[] momentOfInertia(); //diagonal elements of
        // (diagonalized) moment of inertia;
        // should always be a 3-element array
    }

    public interface SphericalTop extends Rotator {
    } //guarantees Ixx = Iyy = Izz

    public interface CylindricalTop extends Rotator {
    } //guarantees Ixx = Iyy

    public interface AsymmetricTop extends Rotator {
    } //all moment-of-inertia elements unequal

}

