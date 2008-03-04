package etomica.atom;

import etomica.api.IAtomType;
import etomica.api.IVector;
import etomica.space.Space;
import etomica.space3d.OrientationFull3D;

/**
 * Molecule class appropriate for a rigid molecule.  The molecule object holds
 * a position and orientation fields.
 *
 * @author Andrew Schultz
 */
public class MoleculeOriented extends Molecule implements IAtomOriented {

    public MoleculeOriented(Space space, IAtomType type) {
        super(type);
        orientation = new OrientationFull3D(space);
        position = space.makeVector();
    }

    public OrientationFull3D getOrientation() {
        return orientation;
    }

    public IVector getPosition() {
        return position;
    }

    private static final long serialVersionUID = 1L;
    protected final OrientationFull3D orientation;
    protected final IVector position;
}
