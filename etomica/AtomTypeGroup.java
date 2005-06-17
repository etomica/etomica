/*
 * Created on May 20, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica;

import etomica.atom.AtomPositionDefinition;
import etomica.utility.Arrays;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class AtomTypeGroup extends AtomType {

    /**
     * used only to create root type
     */
    AtomTypeGroup(AtomIndexManager indexManager) {
        super(indexManager);
    }
    
    /**
     * @param parentType
     * @param positionDefinition
     */
    public AtomTypeGroup(AtomTypeGroup parentType,
            AtomPositionDefinition positionDefinition) {
        super(parentType, positionDefinition);
        if (parentType != null) {
            parentType.childTypes = (AtomType[])Arrays.addObject(parentType.childTypes,this);
        }
    }
    
    /**
     * Returns false, indicating that this is a group type.
     */
    public boolean isLeaf() {
        return false;
    }

    public AtomType[] getChildTypes() {
        return (AtomType[])childTypes.clone();
    }
    
    /**
     * Returns a specified atomType descended from this one in the atomType tree.  
     * Each index of the given array specifies the i-th child at the
     * depth of the array index.  So if path is {2, 0, 3},
     * returns the 3rd childType of the 0th childType of the 2nd childType of
     * this type.  That is: (this type) -> (2nd type) -> (0th type) -> (3rd type)
     */
    public AtomType getDescendant(int[] path) {
        return getDescendant(0, path);
    }
    
    private AtomType getDescendant(int n, int[] path) {
        AtomType child = childTypes[path[n]];
        if(path.length - 1 > n) {//go further down hierarchy
            if(child.isLeaf()) {//no more there
                throw new IllegalArgumentException("Depth of requested descendant exceeds depth of atom hierarchy");
            }//get indicated descendant recursively
            child = ((AtomTypeGroup)child).getDescendant(n+1, path);
        }
        return child;
    }

    
    private AtomType[] childTypes = new AtomType[0];

}
