package etomica.atom;

import etomica.api.IAtom;
import etomica.math.geometry.Shape;


/**
 * Filter that accepts atom if it is inside a specified Shape instance.
 * Position of atom is determined by an AtomPositionDefinition, which
 * if unspecified defaults for each atom to that given by the atom's
 * type. 
 *
 * @author David Kofke
 */

public class AtomFilterInShape implements AtomFilter, java.io.Serializable {

    /**
     * Create filter in which position definition is given by atom's type.
     */
    public AtomFilterInShape(Shape shape) {
        super();
        this.shape = shape;
    }

    /**
     * Returns true if the atom's position is inside the shape.
     */
    public boolean accept(IAtom atom) {
        if(positionDefinition == null) {
            return shape.contains(atom.getType().getPositionDefinition().position(atom));
        }
        return shape.contains(positionDefinition.position(atom));
    }

    /**
     * @return Returns the shape.
     */
    public Shape getShape() {
        return shape;
    }
    /**
     * @param shape The shape to set.
     */
    public void setShape(Shape shape) {
        this.shape = shape;
    }
    /**
     * @return Returns the positionDefinition.
     */
    public AtomPositionDefinition getPositionDefinition() {
        return positionDefinition;
    }
    /**
     * Sets position definition.  If given null, positionDefinition from
     * each atom's type will be used.
     * @param positionDefinition The positionDefinition to set.
     */
    public void setPositionDefinition(AtomPositionDefinition positionDefinition) {
        this.positionDefinition = positionDefinition;
    }
    
    private static final long serialVersionUID = 1L;
    private Shape shape;
    private AtomPositionDefinition positionDefinition;
}
