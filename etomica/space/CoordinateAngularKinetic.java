package etomica.space;

import etomica.Space;


/**
 * Implementation of a coordinate that associates position and orientational
 * coordinates and momenta with each atom.
 */

/*
 * History
 * Created on Jan 26, 2005 by kofke
 */
public class CoordinateAngularKinetic extends CoordinateAngular implements
        ICoordinateAngularKinetic {

    /**
     * Constructs coordinate, making Vector and Orientation instances
     * from the given Space.
     */
    public CoordinateAngularKinetic(Space space) {
        super(space);
        v = space.makeVector();
        omega = space.makeVector();
    }

    /**
     * Returns the angular-velocity vector.
     */
    public Vector angularVelocity() {
        return omega;
    }

    /**
     * Returns the velocity vector.
     */
    public Vector velocity() {
        return v;
    }

    private Vector v, omega;
}
