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
public class Coordinate implements ICoordinate {

    /**
     * 
     */
    public Coordinate(Space space) {
        r = space.makeVector();
    }

    /* (non-Javadoc)
     * @see etomica.space.ICoordinate#position()
     */
    public Vector position() {
        return r;
    }

    protected final Vector r;
}
