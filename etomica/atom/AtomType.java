package etomica.atom;

import etomica.api.IAtomType;
import etomica.api.ISpecies;

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

public abstract class AtomType implements java.io.Serializable, Comparable, IAtomType {

    protected int index;

    protected AtomPositionDefinition positionDefinition;

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

    /* (non-Javadoc)
	 * @see etomica.atom.IAtomType#setIndex(int)
	 */
    public void setIndex(int newIndex) {
        index = newIndex;
    }

    /* (non-Javadoc)
	 * @see etomica.atom.IAtomType#getIndex()
	 */
    public int getIndex() {
        return index;
    }
    
    /* (non-Javadoc)
	 * @see etomica.atom.IAtomType#getSpecies()
	 */
    public abstract ISpecies getSpecies();
    
    /* (non-Javadoc)
	 * @see etomica.atom.IAtomType#getPositionDefinition()
	 */
    public AtomPositionDefinition getPositionDefinition() {
        return positionDefinition;
    }

    /* (non-Javadoc)
	 * @see etomica.atom.IAtomType#setPositionDefinition(etomica.atom.AtomPositionDefinition)
	 */
    public void setPositionDefinition(AtomPositionDefinition newPositionDefinition) {
        positionDefinition = newPositionDefinition;
    }

    /* (non-Javadoc)
	 * @see etomica.atom.IAtomType#setInteracting(boolean)
	 */
    public void setInteracting(boolean b) {
        isInteracting = b;
    }

    /* (non-Javadoc)
	 * @see etomica.atom.IAtomType#isInteracting()
	 */
    public boolean isInteracting() {
        return isInteracting;
    }
    
    /* (non-Javadoc)
	 * @see etomica.atom.IAtomType#compareTo(java.lang.Object)
	 */
    public int compareTo(Object otherAtomType) {
        int otherIndex = ((IAtomType)otherAtomType).getIndex();
        return otherIndex > index ? -1 : (otherIndex == index ? 0 : 1);
    }
    
    /* (non-Javadoc)
	 * @see etomica.atom.IAtomType#toString()
	 */
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

