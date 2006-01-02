package etomica.space;


/**
 * Implemention of a coordinate that has a position and a velocity.
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
     * Set this coordinate's parameters equal to those of the
     * given coordinate.  Overrides superclass to ensure that
     * orientation is copied.  
     * 
     * @throws ClassCastException if argument is not an instance of ICoordinateKinetic
     */
    public void E(ICoordinate coord) {
        this.E((ICoordinateKinetic)coord);
    }
    
    /**
     * Set this coordinate's parameters equal to those of the
     * given coordinate.
     */
    public void E(ICoordinateKinetic coord) {
        r.E(coord.position());
        v.E(coord.velocity());
    }
    
   /**
     * Returns the instance of the velocity vector (not a copy).
     */
    public final Vector velocity() {
        return v;
    }

    protected final Vector v;
    private static final long serialVersionUID = 1L;
}
