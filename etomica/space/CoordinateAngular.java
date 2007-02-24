package etomica.space;

import java.io.Serializable;

/**
 * Implementation of a Coordinate that associates a position and
 * and orientation with an atom, both made by an arbitrary-dimension
 * Space.
 */

public class CoordinateAngular extends Coordinate implements ICoordinateAngular, Serializable {

    /**
     * Makes theh coordinate using Vector and Orientation from the given Space.
     */
    public CoordinateAngular(Space space) {
        super(space);
        orientation = space.makeOrientation();
    }
    
    /**
     * Set this coordinate's parameters equal to those of the
     * given coordinate.  Overrides superclass to ensure that
     * orientation is copied.  
     * 
     * @throws ClassCastException if argument is not an instance of CoordinateAngular
     */
    public void E(ICoordinate coord) {
        super.E(coord);
        orientation.E(((CoordinateAngular)coord).orientation);
    }

    /**
     * Returns the orientation (not a copy).
     */
    public Orientation getOrientation() {
        return orientation;
    }

    protected Orientation orientation;
    private static final long serialVersionUID = 1L;

}
