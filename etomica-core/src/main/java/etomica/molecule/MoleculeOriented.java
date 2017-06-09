/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.api.ISpecies;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.IOrientation3D;
import etomica.space3d.Orientation3D;
import etomica.space3d.OrientationFull3D;

/**
 * Molecule class appropriate for a rigid molecule.  The molecule object holds
 * a position and orientation fields.
 *
 * @author Andrew Schultz
 */
public class MoleculeOriented extends Molecule implements IMoleculeOriented {

    public MoleculeOriented(Space space, ISpecies species, int numLeafAtoms) {
        this(space, species, numLeafAtoms, false);
    }

    public MoleculeOriented(Space space, ISpecies species, int numLeafAtoms, boolean isAxisSymmetric) {
        super(species, numLeafAtoms);
        orientation = isAxisSymmetric ? new Orientation3D(space) : new OrientationFull3D(space);
        position = space.makeVector();
    }

    public IOrientation3D getOrientation() {
        return orientation;
    }

    public Vector getPosition() {
        return position;
    }

    private static final long serialVersionUID = 1L;
    protected final IOrientation3D orientation;
    protected final Vector position;
}
