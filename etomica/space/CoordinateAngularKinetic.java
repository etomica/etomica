package etomica.space;



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
     * Set this coordinate's parameters equal to those of the
     * given coordinate.  Overrides superclass to ensure that
     * velocities are copied.  
     * 
     * @throws ClassCastException if argument is not an instance of CoordinateAngularKinetic
     */
    public void E(ICoordinate coord) {
        super.E(coord);
        v.E(((CoordinateAngularKinetic)coord).v);
        omega.E(((CoordinateAngularKinetic)coord).omega);
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
