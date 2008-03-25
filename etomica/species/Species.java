package etomica.species;

import etomica.api.IAtomPositionDefinition;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IConformation;
import etomica.api.ISpecies;
import etomica.atom.AtomType;
import etomica.simulation.SpeciesManager;
import etomica.util.Arrays;

/**
 * Type for atom that is a group of other atoms, and for which its node is an
 * instance of AtomTreeNodeGroup.
 * 
 * @author andrew
 */
public abstract class Species extends AtomType implements ISpecies {

    /**
     * Simple invokes parent constructor with same arguments.
     */
    public Species(IAtomPositionDefinition positionDefinition) {
        super(positionDefinition);
    }
    
    /**
     * Sets the SpeciesManager.  This is used for callbacks for notification of
     * removal and addition of child types (not that that should ever happen!)
     */
    public void setSpeciesManager(SpeciesManager newSpeciesManager) {
        speciesManager = newSpeciesManager;
        for (int i=0; i<childTypes.length; i++) {
            childTypes[i].setIndex(speciesManager.requestTypeIndex());
        }
    }

    /* (non-Javadoc)
     * @see etomica.atom.IAtomTypeMolecule#removeChildType(etomica.atom.AtomTypeLeaf)
     */
    public void removeChildType(IAtomTypeLeaf removedType) {
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
        childTypes = (IAtomTypeLeaf[])Arrays.removeObject(childTypes,removedType);
        for (int i = 0; i < childTypes.length; i++) {
            childTypes[i].setChildIndex(i);
        }
        if (speciesManager != null) {
            System.err.println("removing child types is generally a bit scary, but you did it while the molecule type was in the simulation, which makes you a bad person");
            speciesManager.atomTypeRemovedNotify(removedType);
        }
    }
    
    /**
     * @return Returns the species.
     */
    public ISpecies getSpecies() {
        return this;
    }

    /* (non-Javadoc)
     * @see etomica.atom.IAtomTypeMolecule#getChildTypes()
     */
    public IAtomTypeLeaf[] getChildTypes() {
        return childTypes;
    }

    /* (non-Javadoc)
     * @see etomica.atom.IAtomTypeMolecule#addChildType(etomica.atom.AtomTypeLeaf)
     */
    public void addChildType(IAtomTypeLeaf newChildType) {
        if (newChildType.getSpecies() != null) {
            throw new IllegalArgumentException(newChildType+" already has a parent");
        }
        newChildType.setSpecies(this);
        newChildType.setChildIndex(childTypes.length);
        childTypes = (IAtomTypeLeaf[]) Arrays.addObject(childTypes, newChildType);
        if (speciesManager != null) {
            System.err.println("You really shouldn't be adding leaf atom types after the Species has been added to the simulation");
            speciesManager.atomTypeAddedNotify(newChildType);
            newChildType.setIndex(speciesManager.requestTypeIndex());
        }
    }
    
    /* (non-Javadoc)
     * @see etomica.atom.IAtomTypeMolecule#setConformation(etomica.config.Conformation)
     */
    public void setConformation(IConformation config) {
        conformation = config;
    }
    
    /* (non-Javadoc)
     * @see etomica.atom.IAtomTypeMolecule#getConformation()
     */
    public IConformation getConformation() {return conformation;}
    
    private static final long serialVersionUID = 2L;
    protected IConformation conformation;
    protected SpeciesManager speciesManager;
    protected IAtomTypeLeaf[] childTypes = new IAtomTypeLeaf[0];
}