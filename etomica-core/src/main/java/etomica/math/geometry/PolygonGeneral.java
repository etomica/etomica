/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.geometry;

import com.fasterxml.jackson.annotation.JsonIgnore;
import etomica.space.Vector;
import etomica.exception.MethodNotImplementedException;
import etomica.space.Space;

/**
 * TODO To change the template for this generated type comment go to Window -
 * Preferences - Java - Code Style - Code Templates
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on May 13, 2005 by kofke
 */
public class PolygonGeneral extends Polygon {

    /**
     * @param embeddedSpace
     * @param nSides
     */
    public PolygonGeneral(Space embeddedSpace, int nSides) {
        super(embeddedSpace, nSides);
    }

    /**
     * @param edges
     */
    public PolygonGeneral(LineSegment[] edges) {
        super(edges);
    }

    /**
     * Returns the value of the area enclosed by the polygon
     */
    //must override in subclass (until a general algorithm is implemented)
    @JsonIgnore
    public double getArea() {
        throw new MethodNotImplementedException(
                "General formula for area not implemented");
    }

    /**
     * Returns true if the given point lies inside or on an edge of the polygon
     */
    //must override in subclass (until a general algorithm is implemented)
    public boolean contains(Vector vector) {
        throw new MethodNotImplementedException(
                "General formula for 'contains' not implemented");
    }

    public void updateVertices() {
        //does nothing, becuse in this case the vertices are the official
        //representation of the polytope
    }

}
