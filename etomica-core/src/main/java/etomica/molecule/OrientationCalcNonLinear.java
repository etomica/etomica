/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.molecule;

import etomica.box.Box;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.IOrientationFull3D;
import etomica.space3d.OrientationFull3D;
import etomica.space3d.RotationTensor3D;

public abstract class OrientationCalcNonLinear implements OrientationCalc {

    @Override
    public int getDOF(Space space) {
        return 3;
    }

    @Override
    public IOrientation makeOrientation(Space space) {
        return new OrientationFull3D(space);
    }

    @Override
    public Vector bodyToSpace(IMolecule molecule, Space space, Vector v) {
        Vector out = space.makeVector();
        out.E(v);
        RotationTensor3D rotationTensor = new RotationTensor3D();
        OrientationFull3D orientation = new OrientationFull3D(space);
        calcOrientation(molecule, orientation);
        rotationTensor.setOrientation(orientation);
        rotationTensor.invert();
        //now rotate back to get space-fixed angular velocity
        rotationTensor.transform(out);
        return out;
    }

    @Override
    public Vector spaceToBody(IMolecule molecule, Space space, Vector v) {
        Vector out = space.makeVector();
        out.E(v);
        RotationTensor3D rotationTensor = new RotationTensor3D();
        OrientationFull3D orientation = new OrientationFull3D(space);
        calcOrientation(molecule, orientation);
        rotationTensor.setOrientation(orientation);
        //now rotate back to get space-fixed angular velocity
        rotationTensor.transform(out);
        return out;
    }

    @Override
    public Vector angularVelocityToMomentum(IMolecule molecule, Box box, Vector omega) {
        Vector L = box.getSpace().makeVector();

        RotationTensor3D rotationTensor = new RotationTensor3D();
        OrientationFull3D orientation = new OrientationFull3D(box.getSpace());
        calcOrientation(molecule, orientation);
        L.E(omega);
        rotationTensor.setOrientation(orientation);
        //find body-fixed angular momentum
        rotationTensor.transform(L);
        //now divide out moment of inertia to get body-fixed angular velocity
        omega.TE(getMomentOfInertia(molecule));
        rotationTensor.invert();
        //now rotate back to get space-fixed angular velocity
        rotationTensor.transform(omega);
        return L;
    }

    @Override
    public Vector angularMomentumToVelocity(IMolecule molecule, Space space, Vector L) {

        OrientationFull3D orientation = new OrientationFull3D(space);
        calcOrientation(molecule, orientation);
        return angularMomentumToVelocity(orientation, molecule, space, L);
    }

    @Override
    public Vector angularMomentumToVelocity(IOrientation orientation, IMolecule molecule, Space space, Vector L) {
        Vector omega = space.makeVector();
        omega.E(L);
        RotationTensor3D rotationTensor = new RotationTensor3D();
        rotationTensor.setOrientation((IOrientationFull3D) orientation);
        //find body-fixed angular momentum
        rotationTensor.transform(omega);
        //now divide out moment of inertia to get body-fixed angular velocity
        omega.DE(molecule.getType().getMomentOfInertia());
        rotationTensor.invert();
        //now rotate back to get space-fixed angular velocity
        rotationTensor.transform(omega);
        return omega;
    }
}
