package etomica.atom;

import etomica.api.IAtomTypeLeaf;
import etomica.api.ISpecies;
import etomica.chem.elements.Element;
import etomica.units.Dimension;
import etomica.units.Mass;

/**
 * Type for an atom that is a leaf in the species hierarchy. An atom of this
 * type is typically representing a physical atom, rather than a group of other
 * atoms.
 * 
 * @author andrew
 */

public class AtomTypeLeaf extends AtomType implements IAtomTypeLeaf {

    public AtomTypeLeaf(Element element) {
        super();
        this.element = element;
        index = -1;
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
    public final Element getElement() {
        return element;
    }

    private static final long serialVersionUID = 1L;
    protected final Element element;
    protected ISpecies parentType;
    protected int childIndex;
}