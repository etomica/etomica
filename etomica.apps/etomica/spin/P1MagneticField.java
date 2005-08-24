package etomica.spin;

import etomica.Space;
import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.potential.Potential1;
import etomica.space.Vector;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on May 24, 2005 by kofke
 */
public class P1MagneticField extends Potential1 {

    /**
     * @param space
     */
    public P1MagneticField(Space space) {
        super(space);
        direction = space.makeVector();
        direction.E(0.0);
        direction.setX(0,1.0);
    }

    /* (non-Javadoc)
     * @see etomica.Potential#energy(etomica.AtomSet)
     */
    public double energy(AtomSet atoms) {
        Vector r = ((Atom)atoms).coord.position();
        return h * r.dot(direction);
    }
    
    
    /**
     * @return Returns the direction.
     */
    public Vector getDirection() {
        return (Vector)direction.clone();
    }
    /**
     * @param direction The direction to set.
     */
    public void setDirection(Vector direction) {
        this.direction.E(direction);
        this.direction.normalize();
    }
    /**
     * @return Returns the h.
     */
    public double getH() {
        return h;
    }
    /**
     * @param h The h to set.
     */
    public void setH(double h) {
        this.h = h;
    }
    private double h;
    private final Vector direction;

}
