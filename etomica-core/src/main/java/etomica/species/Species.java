/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

import etomica.atom.AtomType;
import etomica.config.IConformation;
import etomica.molecule.IMolecule;
import etomica.util.Arrays;

/**
 * Type for atom that is a group of other atoms, and for which its node is an
 * instance of AtomTreeNodeGroup.
 * 
 * @author andrew
 */
public abstract class Species implements ISpecies, java.io.Serializable {

    private static final long serialVersionUID = 2L;
    protected int index;
    protected IConformation conformation;
    protected AtomType[] childTypes = new AtomType[0];

    /**
     * Simple invokes parent constructor with same arguments.
     */
    public Species() {
        super();
    }

    /* (non-Javadoc)
     * @see etomica.atom.IAtomType#getIndex()
     */
    public int getIndex() {
        return index;
    }

    /* (non-Javadoc)
     * @see etomica.atom.IAtomType#setIndex(int)
     */
    public void setIndex(int newIndex) {
        index = newIndex;
    }

    /* (non-Javadoc)
     * @see etomica.atom.IAtomTypeMolecule#removeChildType(etomica.atom.AtomType)
     */
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
        childTypes = (AtomType[]) Arrays.removeObject(childTypes, removedType);
    }

    public AtomType getAtomType(int index) {
        return childTypes[index];
    }

    public int getAtomTypeCount() {
    	return childTypes.length;
    }

    /* (non-Javadoc)
     * @see etomica.atom.IAtomTypeMolecule#addChildType(etomica.atom.AtomType)
     */
    public void addChildType(AtomType newChildType) {
        if (newChildType.getSpecies() != null) {
            throw new IllegalArgumentException(newChildType+" already has a parent");
        }
        newChildType.setSpecies(this);
        childTypes = (AtomType[]) Arrays.addObject(childTypes, newChildType);
    }

    /* (non-Javadoc)
     * @see etomica.atom.IAtomTypeMolecule#getConformation()
     */
    public IConformation getConformation() {
        return conformation;
    }

    /* (non-Javadoc)
     * @see etomica.atom.IAtomTypeMolecule#setConformation(etomica.config.Conformation)
     */
    public void setConformation(IConformation config) {
        conformation = config;
    }

    public void initializeConformation(IMolecule molecule) {
        conformation.initializePositions(molecule.getChildList());
    }
}
