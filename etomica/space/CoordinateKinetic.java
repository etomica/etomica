package etomica.space;

import etomica.Space;

/**
 * Implemention of a coordinate that has a position and a velocity.
  */

/*
 * History
 * Created on Jan 26, 2005 by kofke
 */
public class CoordinateKinetic extends Coordinate implements ICoordinateKinetic {

    /**
     * Constructs object with position and velocity vectors
     * made by the given space.
     */
    public CoordinateKinetic(Space space) {
        super(space);
        v = space.makeVector();
    }
    
    /**
     * Returns the instance of the velocity vector.
     */
    public final Vector velocity() {
        return v;
    }

    protected final Vector v;
}
