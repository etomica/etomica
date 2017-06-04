/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.api.IElement;
import etomica.api.ISpecies;
import etomica.units.Dimension;
import etomica.units.Mass;

/**
 * Type for an atom that is a leaf in the species hierarchy. An atom of this
 * type is typically representing a physical atom, rather than a group of other
 * atoms.
 * 
 * @author andrew
 */

public class AtomTypeLeaf implements IAtomType, Comparable<IAtomType> {

    public AtomTypeLeaf(IElement element) {
        super();
        this.element = element;
        index = -1;
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
     * @see etomica.atom.IAtomTypeLeaf#getParentType()
     */
    public ISpecies getParentType() {
        return parentType;
    }
    
    /* (non-Javadoc)
     * @see etomica.atom.IAtomTypeLeaf#setParentType(etomica.atom.AtomTypeMolecule)
     */
    public void setSpecies(ISpecies newParent) {
        parentType = newParent;
    }
    
    /* (non-Javadoc)
     * @see etomica.atom.IAtomTypeLeaf#setChildIndex(int)
     */
    public void setChildIndex(int newChildIndex) {
        childIndex = newChildIndex;
    }
    
    /* (non-Javadoc)
     * @see etomica.atom.IAtomTypeLeaf#getChildIndex()
     */
    public int getChildIndex() {
        return childIndex;
    }
    
    /* (non-Javadoc)
     * @see etomica.atom.IAtomTypeLeaf#getSpecies()
     */
    public ISpecies getSpecies() {
        return parentType;
    }

    /* (non-Javadoc)
     * @see etomica.atom.IAtomTypeLeaf#getMass()
     */
    public final double getMass() {
        return element.getMass();
    }

    /* (non-Javadoc)
     * @see etomica.atom.IAtomTypeLeaf#rm()
     */
    public final double rm() {
        return element.rm();
    }

    /* (non-Javadoc)
     * @see etomica.atom.IAtomTypeLeaf#getMassDimension()
     */
    public final Dimension getMassDimension() {
        return Mass.DIMENSION;
    }
    
    /* (non-Javadoc)
     * @see etomica.atom.IAtomTypeLeaf#getElement()
     */
    public final IElement getElement() {
        return element;
    }
    
    public int compareTo(IAtomType otherAtomType) {
        int otherIndex = otherAtomType.getIndex();
        return otherIndex > index ? -1 : (otherIndex == index ? 0 : 1);
    }

    private static final long serialVersionUID = 1L;
    protected int index;
    protected final IElement element;
    protected ISpecies parentType;
    protected int childIndex;
}
