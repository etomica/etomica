package etomica.atom;

import etomica.species.ISpecies;

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

    protected int index;

    private AtomPositionDefinition positionDefinition;

    private boolean isInteracting = false;

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
        index = -1;
//        setParentType(null);
    }

    public void setIndex(int newIndex) {
        index = newIndex;
    }

    public int getIndex() {
        return index;
    }
    
    public abstract ISpecies getSpecies();
    
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
    
    public int compareTo(Object otherAtomType) {
        int otherIndex = ((AtomType)otherAtomType).getIndex();
        return otherIndex > index ? -1 : (otherIndex == index ? 0 : 1);
    }
    
    public String toString() {
        return "AtomType "+index;
    }
    
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

