/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.math.geometry.Shape;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculePositionDefinition;


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
        return shape.contains(atom.getPosition());
    }

    public boolean accept(IMolecule mole) {
        return false;
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
    public IMoleculePositionDefinition getPositionDefinition() {
        return positionDefinition;
    }
    /**
     * Sets position definition.  If given null, positionDefinition from
     * each atom's type will be used.
     * @param positionDefinition The positionDefinition to set.
     */
    public void setPositionDefinition(IMoleculePositionDefinition positionDefinition) {
        this.positionDefinition = positionDefinition;
    }
    
    private static final long serialVersionUID = 1L;
    private Shape shape;
    private IMoleculePositionDefinition positionDefinition;
}
