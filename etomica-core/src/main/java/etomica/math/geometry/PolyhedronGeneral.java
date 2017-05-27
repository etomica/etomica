/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.geometry;
import etomica.space.Vector;
import etomica.exception.MethodNotImplementedException;

/**
 * Representation of a mathematical polyhedron, a 3-dimensional polytope.
 */
public class PolyhedronGeneral extends Polyhedron {
    
    protected PolyhedronGeneral(Polygon[] faces) {
        super(faces);
    }
    
    /**
     * Calculates vertices from their internal representation, which
     * for a general polyhedron does nothing because the vertices are
     * the internal representation.
     */
    public void updateVertices() {
        //does nothing, becuse in this case the vertices are the official
        //representation of the polytope
    }
    
    /**
     * Returns the 3-D volume of the polyhedron, which is its conventional volume
     */
    //must override in subclass (until a general algorithm is implemented)
    public double getVolume() {
        throw new MethodNotImplementedException();
    }
    
    /**
     * Returns true if the given point lies inside or on an edge of the polyhedron
     */
    public boolean contains(Vector vector) {
        throw new MethodNotImplementedException("General formula for 'contains' not in place");
    }

}
