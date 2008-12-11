package etomica.species;

import etomica.api.IAtomPositionDefinition;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IConformation;
import etomica.api.ISpecies;
import etomica.atom.AtomType;
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
        super();
        this.positionDefinition = positionDefinition;
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
    }
    
    /**
     * @return Returns the species.
     */
    public ISpecies getSpecies() {
        return this;
    }

    public IAtomTypeLeaf getChildType(int index) {
    	return childTypes[index];
    }

    public int getChildTypeCount() {
    	return childTypes.length;
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
 
    
    /* (non-Javadoc)
     * @see etomica.atom.IAtomType#getPositionDefinition()
     */
    public IAtomPositionDefinition getPositionDefinition() {
        return positionDefinition;
    }

    /* (non-Javadoc)
     * @see etomica.atom.IAtomType#setPositionDefinition(etomica.atom.AtomPositionDefinition)
     */
    public void setPositionDefinition(IAtomPositionDefinition newPositionDefinition) {
        positionDefinition = newPositionDefinition;
    }

    private static final long serialVersionUID = 2L;
    protected IConformation conformation;
    protected IAtomPositionDefinition positionDefinition;
    protected IAtomTypeLeaf[] childTypes = new IAtomTypeLeaf[0];
}