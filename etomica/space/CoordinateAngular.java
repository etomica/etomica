package etomica.space;



/**
 * Implementation of a coordinate that associates a position and
 * and orientation with an atom, both made by an arbitrary-dimension
 * Space.
  */

/*
 * History
 * Created on Jan 26, 2005 by kofke
 */
public class CoordinateAngular extends Coordinate implements ICoordinateAngular {

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

    /* (non-Javadoc)
     * @see etomica.space.ICoordinateAngular#orientation()
     */
    public Orientation orientation() {
        return orientation;
    }

    protected Orientation orientation;
}
