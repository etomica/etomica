package etomica.atom;

import etomica.api.ISpecies;
import etomica.api.IVectorMutable;
import etomica.space.ISpace;
import etomica.space3d.OrientationFull3D;

/**
 * Molecule class appropriate for a rigid molecule.  The molecule object holds
 * a position and orientation fields.
 *
 * @author Andrew Schultz
 */
public class MoleculeOriented extends Molecule implements IAtomOriented {

    public MoleculeOriented(ISpace space, ISpecies species, int numLeafAtoms) {
        super(species, numLeafAtoms);
        orientation = new OrientationFull3D(space);
        position = space.makeVector();
    }

    public OrientationFull3D getOrientation() {
        return orientation;
    }

    public IVectorMutable getPosition() {
        return position;
    }

    private static final long serialVersionUID = 1L;
    protected final OrientationFull3D orientation;
    protected final IVectorMutable position;
}
