package etomica.atom;

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

public class AtomTypeLeaf extends AtomType {

    public AtomTypeLeaf(Element element) {
        this(element, new AtomPositionDefinitionSimple());
    }
    
    /**
     * Invokes parent constructor with and sets the mass to the given value.
     */
    public AtomTypeLeaf(Element element,
            AtomPositionDefinition positionDefinition) {
        super(positionDefinition);
        this.element = element;
        index = -1;
    }

    /**
     * Returns the parent AtomType of this AtomType.
     */
    public AtomTypeMolecule getParentType() {
        return parentType;
    }
    
    public void setParentType(AtomTypeMolecule newParent) {
        parentType = newParent;
    }
    
    public void setChildIndex(int newChildIndex) {
        childIndex = newChildIndex;
    }
    
    public int getChildIndex() {
        return childIndex;
    }
    
    public ISpecies getSpecies() {
        return parentType.getSpecies();
    }

    /**
     * Returns the value of the mass.
     */
    public final double getMass() {
        return element.getMass();
    }

    /**
     * Returns the reciprocal of the mass, 1.0/mass
     */
    public final double rm() {
        return element.rm();
    }

    /**
     * Returns Dimension.MASS, indicating that "mass" has dimensions of mass.
     */
    public final Dimension getMassDimension() {
        return Mass.DIMENSION;
    }
    
    public final Element getElement() {
        return element;
    }

    protected final Element element;
    protected AtomTypeMolecule parentType;
    protected int childIndex;
}