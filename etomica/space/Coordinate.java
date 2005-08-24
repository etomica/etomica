package etomica.space;



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
public class Coordinate implements ICoordinate, java.io.Serializable {

    /**
     * Makes the coordinate vector using the given Space.
     */
    public Coordinate(Space space) {
        r = space.makeVector();
    }
    
    
    /**
     * Set this coordinate's parameters equal to those of the
     * given coordinate.
     */
    public void E(ICoordinate coord) {
        r.E(coord.position());
    }

    /**
     * Returns the position vector (not a copy).
     */
    public final Vector position() {
        return r;
    }

    protected final Vector r;
}
