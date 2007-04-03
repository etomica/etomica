package etomica.atom;

import etomica.species.Species;
import etomica.util.Arrays;

/**
 * Type for atom that is a group of other atoms, and for which its node is an
 * instance of AtomTreeNodeGroup.
 * 
 * @author andrew
 */
public class AtomTypeGroup extends AtomType {

    /**
     * Used only to create root type
     */
    AtomTypeGroup(AtomAddressManager indexManager) {
        super(indexManager);
        assignChildOrdinals();
    }

    /**
     * Simple invokes parent constructor with same arguments.
     */
    public AtomTypeGroup(AtomPositionDefinition positionDefinition) {
        super(positionDefinition);
        assignChildOrdinals();
    }

    public void resetIndex() {
        boolean recurse = index == -1;
        super.resetIndex();
        if (recurse) {
            for (int i=0; i<childTypes.length; i++) {
                childTypes[i].resetIndex();
            }
        }
    }

    public void removeChildType(AtomType removedType) {
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
        childTypes = (AtomType[])Arrays.removeObject(((AtomTypeGroup)childTypes[0]).childTypes,removedType);
        assignChildOrdinals();
        if (parentType != null) {
            parentType.childTypeRemovedNotify(removedType);
        }
    }
    
    void childTypeRemovedNotify(AtomType removedType) {
        if (parentType != null) {
            parentType.childTypeRemovedNotify(removedType);
        }
    }

    /**
     * Assigns ordinals to all child atoms, numbering them sequentially
     * according to their position in the childList.
     */
    private void assignChildOrdinals() {
        for (int i = 0; i < childTypes.length; i++) {
            childTypes[i].resetIndex();
        }
    }

    int requestIndex() {
        if(parentType == null) return -1;
        return parentType.requestIndex();
    }
    
    /**
     * Returns false, indicating that this is a group type.
     */
    public boolean isLeaf() {
        return false;
    }

    /**
     * Returns (a copy of) the array of child types of this group.
     */
    public AtomType[] getChildTypes() {
        return (AtomType[]) childTypes.clone();
    }

    /**
     * @param species
     *            The species to set.
     */
    public void setSpecies(Species species) {
        super.setSpecies(species);
        for (int i=0; i<childTypes.length; i++) {
            childTypes[i].setSpecies(species);
        }
    }

    /**
     * Returns a specified atomType descended from this one in the atomType
     * tree. Each index of the given array specifies the i-th child at the depth
     * of the array index. So if path is {2, 0, 3}, returns the 3rd childType of
     * the 0th childType of the 2nd childType of this type. That is: (this type) ->
     * (2nd type) -> (0th type) -> (3rd type)
     */
    public AtomType getDescendant(int[] path) {
        return getDescendant(0, path);
    }

    private AtomType getDescendant(int n, int[] path) {
        AtomType child = childTypes[path[n]];
        if (path.length - 1 > n) {//go further down hierarchy
            if (child.isLeaf()) {//no more there
                throw new IllegalArgumentException(
                        "Depth of requested descendant exceeds depth of atom hierarchy");
            }//get indicated descendant recursively
            child = ((AtomTypeGroup) child).getDescendant(n + 1, path);
        }
        return child;
    }

    void addChildType(AtomType newChildType) {
        if (species != null) {
            newChildType.setSpecies(species);
        }
        childTypes = (AtomType[]) Arrays.addObject(childTypes, newChildType);
        childTypeAddedNotify(newChildType);
    }
    
    protected void childTypeAddedNotify(AtomType newChildType) {
        if (parentType != null) {
            parentType.childTypeAddedNotify(newChildType);
        }
    }

    private static final long serialVersionUID = 1L;
    protected AtomType[] childTypes = new AtomType[0];
}