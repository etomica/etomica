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
public class CoordinateKinetic extends Coordinate implements ICoordinateKinetic {

    /**
     * @param a
     */
    public CoordinateKinetic(Space space) {
        super(space);
        v = space.makeVector();
    }
    
    public Vector velocity() {
        return v;
    }

    protected final Vector v;
}
