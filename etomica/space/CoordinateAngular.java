package etomica.space;

import etomica.Space;


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

    /**
     * @param space
     */
    public CoordinateAngular(Space space) {
        super(space);
        orientation = space.makeOrientation();
    }

    /* (non-Javadoc)
     * @see etomica.space.ICoordinateAngular#orientation()
     */
    public Orientation orientation() {
        return orientation;
    }

    protected Orientation orientation;
}
