package etomica.atom;

import etomica.Parameter;
import etomica.Parameter.Source;
import etomica.species.Species;
import etomica.util.Arrays;
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

    public static Parameter.Source[] parameterSource = new Parameter.Source[0];
    AtomFactory creator;//set in constructor of AtomFactory
    public Parameter[] parameter;
    protected int speciesIndex = -1;
    private Species species;
    private final int index;

    //fields for linked list of all instances of AtomType
    public final AtomType previousInstance;
    private static AtomType lastInstance;

    private final AtomIndexManager indexManager;
    private AtomPositionDefinition positionDefinition;

    private AtomTypeGroup parentType;
    
    private boolean isInteracting = false;

    //    private Parameter.Electrostatic electroParameter;

    /**
     * Used only to create root type.
     */
    AtomType(AtomIndexManager indexManager) {
        parentType = null;
        this.indexManager = indexManager;
        positionDefinition = null;
        previousInstance = null;
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
    public AtomType(AtomTypeGroup parentType,
            AtomPositionDefinition positionDefinition) {
        this.parentType = parentType;
        if (parentType == null) {
            indexManager = AtomIndexManager
                    .makeSimpleIndexManager(Default.BIT_LENGTH);
            index = 0;
        } else {
            indexManager = parentType.getIndexManager().makeChildManager();
            parentType.childTypes = (AtomType[]) Arrays.addObject(
                    parentType.childTypes, this);
            index = parentType.requestIndex();
        }
        this.positionDefinition = positionDefinition;

        //update linked list of instances
        this.previousInstance = lastInstance;
        lastInstance = this;

        //set up global parameters
        parameter = new Parameter[parameterSource.length];
        for (int i = 0; i < parameter.length; i++) {
            parameter[i] = parameterSource[i].makeParameter();
        }
    }
    
    public int getIndex() {
        return index;
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
    public AtomIndexManager getIndexManager() {
        return indexManager;
    }

    /**
     * Implementation of Comparable interface, based on type index value.  
     * Returns -1, 0, 1 if getIndexManager().getTypeIndex() for this type
     * is less, equal, or greater, respectively, than the corresponding 
     * value for the given type.
     */
    public int compareTo(Object atomType) {
        int otherIndex = ((AtomType) atomType).indexManager.getTypeIndex();
        int myIndex = indexManager.getTypeIndex();
        return otherIndex > myIndex ? 1 : (otherIndex == myIndex ? 0 : -1);
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
     * Returns true if an atom of this type is descended from an atom (any atom)
     * having the given type.
     */
    public boolean isDescendedFrom(AtomType type) {
        return indexManager.isDescendedFrom(type.indexManager);
    }

    /**
     * Method not yet implemented. Under developement
     */
    //TODO implement or remove the Parameter facility in AtomType
    protected void addGlobalParameter(Parameter.Source source) {
        Parameter[] newParameter = new Parameter[parameter.length + 1];
        for (int i = 0; i < parameter.length; i++)
            newParameter[i] = parameter[i];
        newParameter[parameter.length] = source.makeParameter();
        parameter = newParameter;
    }

    /**
     * Adds given parameter source to parameter-source array and returns index
     * indicating where in atomtype parameter-array the source's parameter will
     * be placed.
     */
    public static int requestParameterIndex(Parameter.Source source) {
        Parameter.Source[] newSource = new Parameter.Source[parameterSource.length + 1];
        for (int i = 0; i < parameterSource.length; i++)
            newSource[i] = parameterSource[i];
        int index = parameterSource.length;
        newSource[index] = source;
        parameterSource = newSource;

        //make parameter for any existing AtomType instances
        for (AtomType t = lastInstance; t != null; t = t.previousInstance) {
            t.addGlobalParameter(source);
        }
        return index;
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
    public int getDepth() {
        return indexManager.getDepth();
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

