/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;

/**
 * Molecule class appropriate for a rigid molecule in a dynamic context.  The
 * molecule object holds a position, velocity orientation and angular momentum
 * as fields.
 *
 * @author Andrew Schultz
 */
public class MoleculeOrientedDynamic extends MoleculeOriented implements IMoleculeOrientedKinetic {

    public MoleculeOrientedDynamic(Space space, ISpecies species, int numLeafAtoms) {
        super(space, species, numLeafAtoms);
        angularMomentum = space.makeVector();
        velocity = space.makeVector();
    }

    public Vector getAngularVelocity() {
        return angularMomentum;
    }

    public Vector getVelocity() {
        return velocity;
    }

    private static final long serialVersionUID = 1L;
    protected final Vector angularMomentum;
    protected final Vector velocity;
}
