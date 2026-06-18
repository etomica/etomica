/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.molecule;

import etomica.box.Box;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Orientation3D;
import etomica.space3d.RotationTensor3D;

public abstract class OrientationCalcLinear implements OrientationCalc {

    @Override
    public int getDOF(Space space) {
        return space.D() - 1;
    }

    @Override
    public IOrientation makeOrientation(Space space) {
        return space.D() == 3 ? new Orientation3D(space) : space.makeOrientation();
    }

    @Override
    public Vector bodyToSpace(IMolecule molecule, Space space, Vector v) {

        IMolecule molecule0 = molecule.getType().makeMolecule();
        molecule0.getType().getConformation().initializePositions(molecule0.getChildList());
        IOrientation orientation0 = makeOrientation(space);
        calcOrientation(molecule0, orientation0);
        Vector direction0 = orientation0.getDirection();

        IOrientation orientation = makeOrientation(space);
        calcOrientation(molecule, orientation);
        Vector direction = orientation.getDirection();

        RotationTensor3D tensor = new RotationTensor3D();
        tensor.setRotation(direction0, direction);

        Vector rv = space.makeVector();
        rv.E(v);
        tensor.transform(rv);
        return rv;
    }

    @Override
    public Vector spaceToBody(IMolecule molecule, Space space, Vector v) {

        IMolecule molecule0 = molecule.getType().makeMolecule();
        molecule0.getType().getConformation().initializePositions(molecule0.getChildList());
        IOrientation orientation0 = makeOrientation(space);
        calcOrientation(molecule0, orientation0);
        Vector direction0 = orientation0.getDirection();

        IOrientation orientation = makeOrientation(space);
        calcOrientation(molecule, orientation);
        Vector direction = orientation.getDirection();

        RotationTensor3D tensor = new RotationTensor3D();
        tensor.setRotation(direction, direction0);

        Vector rv = space.makeVector();
        rv.E(v);
        tensor.transform(rv);
        return rv;
    }

    @Override
    public Vector angularVelocityToMomentum(IMolecule molecule, Box box, Vector omega) {
        Vector L = box.getSpace().makeVector();
        L.E(omega);
        L.TE(Math.sqrt(0.5*getMomentOfInertia(molecule).squared()));
        return L;
    }

    @Override
    public Vector angularMomentumToVelocity(IMolecule molecule, Space space, Vector L) {
        Vector omega = space.makeVector();
        omega.E(L);
        omega.TE(1.0/Math.sqrt(0.5*getMomentOfInertia(molecule).squared()));
        return omega;
    }

    @Override
    public Vector angularMomentumToVelocity(IOrientation orientation, IMolecule molecule, Space space, Vector L) {
        return angularMomentumToVelocity(molecule, space, L);
    }
}
