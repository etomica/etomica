package etomica.atom;

import etomica.space.IVector;
import etomica.space.Space;
import etomica.space3d.OrientationFull3D;

/**
 * Molecule class appropriate for a rigid molecule in a dynamic context.  The
 * molecule object holds a position, velocity orientation and angular momentum
 * as fields.
 *
 * @author Andrew Schultz
 */
public class MoleculeOrientedDynamic extends MoleculeOriented implements IAtomOrientedKinetic {

    public MoleculeOrientedDynamic(Space space, AtomType type) {
        super(space, type);
        angularMomentum = space.makeVector();
        velocity = space.makeVector();
    }

    public OrientationFull3D getOrientation() {
        return orientation;
    }

    public IVector getAngularVelocity() {
        return angularMomentum;
    }

    public IVector getPosition() {
        return position;
    }

    public IVector getVelocity() {
        return velocity;
    }

    private static final long serialVersionUID = 1L;
    protected final IVector angularMomentum;
    protected final IVector velocity;
}
