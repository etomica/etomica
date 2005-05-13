package etomica.math.geometry;

import etomica.Space;
import etomica.space.Vector;

/**
 * TODO To change the template for this generated type comment go to Window -
 * Preferences - Java - Code Style - Code Templates
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