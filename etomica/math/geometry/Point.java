package etomica.math.geometry;

import etomica.space.Space;
import etomica.space.Vector;

/**
 * A zero-dimensional polytope, a mathematical point in space.
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on May 12, 2005 by kofke
 */
public class Point extends Polytope {

    public Point(Space embeddedSpace) {
        this(embeddedSpace, embeddedSpace.makeVector());
    }

    public Point(Space embeddedSpace, Vector vector) {
        super(embeddedSpace, vector);
    }

    /**
     * Does nothing; the point is the vertex.
     */
    public void updateVertices() {
        //does nothing because vertex is the internal representation of point
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

}