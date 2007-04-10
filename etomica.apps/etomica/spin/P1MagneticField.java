package etomica.spin;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomSet;
import etomica.potential.Potential1;
import etomica.space.IVector;
import etomica.space.Space;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
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
        IVector r = ((AtomLeaf)atoms).getPosition();
        return h * r.dot(direction);
    }
    
    
    /**
     * @return Returns the direction.
     */
    public IVector getDirection() {
        return direction;
    }
    /**
     * @param direction The direction to set.
     */
    public void setDirection(IVector direction) {
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

    private static final long serialVersionUID = 1L;
    private double h;
    private final IVector direction;
}
