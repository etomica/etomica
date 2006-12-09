package etomica.action;
import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.space.ICoordinateKinetic;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * 
 * Changes the velocity of an atom by a specified vector amount.
 * To accelerate all atoms in a molecule (or atom group), wrap an
 * instance of this class in an AtomGroupAction.
 * 
 * @author David Kofke
 */
public class AtomActionAccelerateBy extends AtomActionAdapter {
    
    private static final long serialVersionUID = 1L;
    private final Vector accelerationVector;
    
    public AtomActionAccelerateBy(Space space) {
        accelerationVector = space.makeVector();
    }
    
    public void actionPerformed(Atom atom) {
        ((ICoordinateKinetic)((AtomLeaf)atom).getCoord()).velocity().PE(accelerationVector);
    }
       
    /**
     * Returns the acceleration vector, the acceleration that will be
     * applied by this action. Returns the vector used by this
     * instance, not a copy, so any manipulation of the returned vector will
     * affect the action of this instance.
     */
    public Vector getAccelerationVector() {
        return accelerationVector;
    }
    /**
     * @param accelerationVector The acceleration vector to set.  A local copy
     * is made of the given vector.
     */
    public void setAccelerationVector(Vector accelerationVector) {
        this.accelerationVector.E(accelerationVector);
    }
}