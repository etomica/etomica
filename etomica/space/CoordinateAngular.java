package etomica.space;

import etomica.Space;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
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
