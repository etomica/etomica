package etomica.atom;

import etomica.config.Conformation;
import etomica.simulation.SpeciesManager;
import etomica.species.ISpecies;
import etomica.util.Arrays;

/**
 * Type for atom that is a group of other atoms, and for which its node is an
 * instance of AtomTreeNodeGroup.
 * 
 * @author andrew
 */
public class AtomTypeMolecule extends AtomType {

    /**
     * Simple invokes parent constructor with same arguments.
     */
    public AtomTypeMolecule(AtomPositionDefinition positionDefinition) {
        super(positionDefinition);
    }
    
    /**
     * Sets the Species.  This may only be called once (typically from the
     * Species constructor)
     */
    public void setSpecies(ISpecies species) {
        if (this.species != null) {
            throw new RuntimeException("Can only call setSpecies once");
        }
        this.species = species;
    }
    
    public void setSpeciesManager(SpeciesManager newSpeciesManager) {
        speciesManager = newSpeciesManager;
        index = speciesManager.requestTypeIndex();
        for (int i=0; i<childTypes.length; i++) {
            childTypes[i].setIndex(speciesManager.requestTypeIndex());
        }
    }

    public void removeChildType(AtomTypeLeaf removedType) {
        boolean success = false;
        for (int i=0; i<childTypes.length; i++) {
            if (childTypes[i] == removedType) {
                success = true;
                break;
            }
        }
        if (!success) {
            throw new IllegalArgumentException("AtomType "+removedType+" is not my child!");
        }
        childTypes = (AtomTypeLeaf[])Arrays.removeObject(childTypes,removedType);
        for (int i = 0; i < childTypes.length; i++) {
            childTypes[i].setChildIndex(i);
        }
        if (speciesManager != null) {
            speciesManager.atomTypeRemovedNotify(removedType);
        }
    }
    
    /**
     * @return Returns the species.
     */
    public ISpecies getSpecies() {
        return species;
    }

    /**
     * Returns the array of child types of this group.
     */
    public AtomTypeLeaf[] getChildTypes() {
        return childTypes;
    }

    public void addChildType(AtomTypeLeaf newChildType) {
        if (newChildType.getParentType() != null) {
            throw new IllegalArgumentException(newChildType+" already has a parent");
        }
        newChildType.setParentType(this);
        newChildType.setChildIndex(childTypes.length);
        childTypes = (AtomTypeLeaf[]) Arrays.addObject(childTypes, newChildType);
        if (speciesManager != null) {
            speciesManager.atomTypeAddedNotify(newChildType);
            newChildType.setIndex(speciesManager.requestTypeIndex());
        }
    }
    
    /**
     * Sets the conformation used to set the standard arrangement of
     * the atoms/atom-groups produced by this factory.
     */
    public void setConformation(Conformation config) {
        conformation = config;
    }
    
    /**
     * Returns the conformation used to set the standard arrangement of
     * the atoms/atom-groups produced by this factory.
     */
    public Conformation getConformation() {return conformation;}
    
    private static final long serialVersionUID = 2L;
    protected Conformation conformation;
    protected ISpecies species;
    protected SpeciesManager speciesManager;
    protected AtomTypeLeaf[] childTypes = new AtomTypeLeaf[0];
}