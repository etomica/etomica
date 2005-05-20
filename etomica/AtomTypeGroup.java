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

    public AtomType[] getChildTypes() {
        return (AtomType[])childTypes.clone();
    }
    
    private AtomType[] childTypes = new AtomType[0];

}
