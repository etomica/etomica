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
public class CoordinateAngularKinetic extends CoordinateAngular implements
        ICoordinateAngularKinetic {

    /**
     * @param space
     */
    public CoordinateAngularKinetic(Space space) {
        super(space);
        v = space.makeVector();
        omega = space.makeVector();
    }

    /* (non-Javadoc)
     * @see etomica.space.ICoordinateAngularKinetic#angularVelocity()
     */
    public Vector angularVelocity() {
        return omega;
    }

    public Vector velocity() {
        return v;
    }

    private Vector v, omega;
}
