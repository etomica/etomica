package etomica.action;
import etomica.Atom;
import etomica.Space;
import etomica.data.DataSourceVelocityAverage;
import etomica.space.Vector;

/**
 * Moves (translates) an atom to a specified position.  Location of the
 * atom (which may be an atom group) is determined by an AtomPositionDefinition
 * instance that may be set for this class.
 */
public class AtomActionAccelerateTo extends AtomActionAdapter {
    
    private final Vector dr;
    private final Vector targetVelocity;
    private final AtomGroupAction atomAccelerator;
    private final AtomGroupAction velocityMeter;

    /**
     * Creates new action with atom position defined by its
     * center of mass (via DataSourceCOM).
     * @param space
     */
    public AtomActionAccelerateTo(Space space) {
        dr = space.makeVector();
        targetVelocity = space.makeVector();
        atomAccelerator = new AtomGroupAction(new AtomActionAccelerateBy(space));
        velocityMeter = new AtomGroupAction(new DataSourceVelocityAverage(space));
    }
    
    public void actionPerformed(Atom atom) {
        velocityMeter.actionPerformed(atom);
        Vector currentVelocity = ((DataSourceVelocityAverage)velocityMeter.getAction()).getVelocityAverage();
        dr.Ev1Mv2(targetVelocity, currentVelocity);
        ((AtomActionAccelerateBy)atomAccelerator.getAction()).setAccelerationVector(dr);
        atomAccelerator.actionPerformed(atom);
    }
       
    /**
     * @return Returns the targetVelocity, the velocity that the
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