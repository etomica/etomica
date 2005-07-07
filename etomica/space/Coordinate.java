package etomica.space;

import etomica.Space;


/**
 * Implementation of ICoordinate interface in which
 * the only atom state parameter is its position, which
 * is represented by a Vector from an arbitrary-dimension
 * Space.
 */

/*
 * History
 * Created on Jan 26, 2005 by kofke
 */
public class Coordinate implements ICoordinate {

    /**
     * Makes the coordinate vector using the given Space.
     */
    public Coordinate(Space space) {
        r = space.makeVector();
    }

    /**
     * Returns the position vector (not a copy).
     */
    public final Vector position() {
        return r;
    }

    protected final Vector r;
}
