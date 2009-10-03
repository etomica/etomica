package etomica.atom;

import etomica.api.ISpecies;
import etomica.api.IVectorMutable;
import etomica.space.ISpace;

/**
 * Molecule class appropriate for a rigid molecule in a dynamic context.  The
 * molecule object holds a position, velocity orientation and angular momentum
 * as fields.
 *
 * @author Andrew Schultz
 */
public class MoleculeOrientedDynamic extends MoleculeOriented implements IMoleculeOrientedKinetic {

    public MoleculeOrientedDynamic(ISpace space, ISpecies species, int numLeafAtoms) {
        super(space, species, numLeafAtoms);
        angularMomentum = space.makeVector();
        velocity = space.makeVector();
    }

    public IVectorMutable getAngularVelocity() {
        return angularMomentum;
    }

    public IVectorMutable getVelocity() {
        return velocity;
    }

    private static final long serialVersionUID = 1L;
    protected final IVectorMutable angularMomentum;
    protected final IVectorMutable velocity;
}
