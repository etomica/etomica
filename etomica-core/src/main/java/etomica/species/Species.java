/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

import etomica.atom.IAtomType;
import etomica.api.IMolecule;
import etomica.api.ISpecies;
import etomica.config.IConformation;
import etomica.util.Arrays;

/**
 * Type for atom that is a group of other atoms, and for which its node is an
 * instance of AtomTreeNodeGroup.
 * 
 * @author andrew
 */
public abstract class Species implements ISpecies, java.io.Serializable {

    protected int index;

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

    /**
     * Simple invokes parent constructor with same arguments.
     */
    public Species() {
        super();
    }

    /* (non-Javadoc)
     * @see etomica.atom.IAtomTypeMolecule#removeChildType(etomica.atom.AtomTypeLeaf)
     */
    public void removeChildType(IAtomType removedType) {
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
        childTypes = (IAtomType[])Arrays.removeObject(childTypes,removedType);
        for (int i = 0; i < childTypes.length; i++) {
            childTypes[i].setChildIndex(i);
        }
    }

    public IAtomType getAtomType(int index) {
    	return childTypes[index];
    }

    public int getAtomTypeCount() {
    	return childTypes.length;
    }

    /* (non-Javadoc)
     * @see etomica.atom.IAtomTypeMolecule#addChildType(etomica.atom.AtomTypeLeaf)
     */
    public void addChildType(IAtomType newChildType) {
        if (newChildType.getSpecies() != null) {
            throw new IllegalArgumentException(newChildType+" already has a parent");
        }
        newChildType.setSpecies(this);
        newChildType.setChildIndex(childTypes.length);
        childTypes = (IAtomType[]) Arrays.addObject(childTypes, newChildType);
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
 
    public void initializeConformation(IMolecule molecule) {
        conformation.initializePositions(molecule.getChildList());
    }
    private static final long serialVersionUID = 2L;
    protected IConformation conformation;
    protected IAtomType[] childTypes = new IAtomType[0];
}
