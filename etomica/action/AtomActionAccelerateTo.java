package etomica.action;
import etomica.atom.Atom;
import etomica.atom.AtomGroupVelocityAverage;
import etomica.atom.AtomLeaf;
import etomica.space.ICoordinateKinetic;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Sets the velocity of an atom to a specified vector value.  If applied
 * to a molecule, works with the mass-averaged velocity of the atoms it comprises,
 * and adds a constant velocity to all atoms to make the mass-averaged velocity
 * equal a target value.
 */
public class AtomActionAccelerateTo extends AtomActionAdapter {
    
    private static final long serialVersionUID = 1L;
    private final Vector dr;
    private final Vector targetVelocity;
    private final AtomGroupAction atomAccelerator;
    private final AtomGroupVelocityAverage velocityMeter;

    /**
     * Creates new action with target velocity equal to zero.
     */
    public AtomActionAccelerateTo(Space space) {
        dr = space.makeVector();
        targetVelocity = space.makeVector();
        atomAccelerator = new AtomGroupAction(new AtomActionAccelerateBy(space));
        velocityMeter = new AtomGroupVelocityAverage(space);
    }
    
    /**
     * Adds velocity vector to the atom (or to all atoms forming it) 
     * so that its mass-average velocity equals the target velocity.
     */
    public void actionPerformed(Atom atom) {
        if(atom.getType().isLeaf()) {
            ((ICoordinateKinetic)((AtomLeaf)atom).getCoord()).getVelocity().E(targetVelocity);
        } else {
            Vector currentVelocity = velocityMeter.getVelocityAverage(atom);
            dr.Ev1Mv2(targetVelocity, currentVelocity);
            ((AtomActionAccelerateBy)atomAccelerator.getAction()).setAccelerationVector(dr);
            atomAccelerator.actionPerformed(atom);
        }
    }
       
    /**
     * @return Returns the targetVelocity, the mass-average velocity that the
     * atom will be accelerated to by this action.
     */
    public Vector getTargetVelocity() {
        return targetVelocity;
    }
    /**
     * @param destination The velocity to set.  A local copy
     * is made of the given vector.
     */
    public void setTargetVelocity(Vector targetVelocity) {
        this.targetVelocity.E(targetVelocity);
    }
}
