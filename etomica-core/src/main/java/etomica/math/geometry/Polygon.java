/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.geometry;

import etomica.space.Vector;
import etomica.space.Space;

/**
 * Representation of a mathematical polygon, a 2-dimensional polytope. Contains
 * all information needed to represent a polygon, methods to set its position
 * and orientation, and methods that return an external representation of it. 
 * Provides no external means to change the polygon size and shape, as the form
 * and capabilities of mutator methods are particular to the type of polygon
 * defined by the subclass.
 */
public abstract class Polygon extends Polytope {

    /**
     * Constructs a polygon with the given number of sides arranged in a closed loop.
     */
    protected Polygon(Space embeddedSpace, int nSides) {
        this(makeEdges(embeddedSpace, nSides));
    }

    /**
     * Constructs a polygon using the given edges for its sides.
     */
    protected Polygon(LineSegment[] edges) {
        super(edges);
        this.edges = edges;
    }
    
    private static LineSegment[] makeEdges(Space embeddedSpace, int nSides) {
        Vector[] vertices = embeddedSpace.makeVectorArray(nSides);
        LineSegment[] edges = new LineSegment[nSides];
        for (int i = 1; i < nSides; i++) {
            edges[i] = new LineSegment(embeddedSpace,
                    vertices[i-1], vertices[i]);
        }
        edges[0] = new LineSegment(embeddedSpace, vertices[nSides-1],
                vertices[0]);
        return edges;
    }

    /**
     * Returns all edges defined by the polygon.
     */
    public LineSegment[] getEdges() {
        updateVertices();
        return edges;
    }
    
    /**
     * Returns the 2-D volume of the polygon, which is its area
     */
    public double getVolume() {
        return getArea();
    }
    
    /**
     * Returns the sum of the length of the edges of the polygon
     */
    public double getPerimeter() {
        updateVertices();
        double sum = 0.0;
        for(int i=0; i<edges.length; i++) {
            sum += edges[i].getLength();
        }
        return sum;
    }

    public abstract double getArea();
    
    protected final LineSegment[] edges;

}
