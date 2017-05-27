/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.geometry;

import etomica.space.Vector;
import etomica.space.Space;

/**
 * A zero-dimensional polytope, a mathematical point in space.
 * 
 * @author David Kofke
 *  
 */
public class Point extends Polytope {

    public Point(Space embeddedSpace) {
        this(embeddedSpace, embeddedSpace.makeVector());
    }

    public Point(Space embeddedSpace, Vector vectorRandom) {
        super(embeddedSpace, vectorRandom);
    }

    /**
     * Does nothing; the point is the vertex.
     */
    public void updateVertices() {
        //does nothing because vertex is the internal representation of point
    }
    
    public LineSegment[] getEdges() {
    	return edges;
    }

    /**
     * Returns zero.
     */
    public double getVolume() {
        return 0;
    }

    /**
     * Returns true if this point is at exactly the same location as the point
     * defined by the given vector.
     */
    public boolean contains(Vector v) {
        return vertices[0].equals(v);
    }
    
    public String toString() {
        return vertices[0].toString();
    }

    private static final long serialVersionUID = 1L;
    private final static LineSegment[] edges = new LineSegment[0];
}
