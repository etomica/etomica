/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.box.Box;
import etomica.space.Vector;
import etomica.space3d.IOrientationFull3D;
import etomica.space3d.OrientationFull3D;
import etomica.space3d.RotationTensor3D;

/**
 * OrientationCalc implementation that handles a monotomic oriented molecule.
 *
 * @author Andrew Schultz
 */
public class OrientationCalcAtom implements OrientationCalc {

    public void calcOrientation(IMolecule molecule,
                                IOrientationFull3D orientation) {
        orientation.E(((IMoleculeOriented) molecule).getOrientation());
    }

    public void setOrientation(IMolecule molecule, Box box,
                               IOrientationFull3D orientation) {
        ((IMoleculeOriented) molecule).getOrientation().E(orientation);
    }

    @Override
    public Vector getAngularMomentum(IMolecule molecule, Vector com, Box box) {
        Vector omega = ((IMoleculeOrientedKinetic) molecule).getAngularVelocity();
        Vector L = box.getSpace().makeVector();

        RotationTensor3D rotationTensor = new RotationTensor3D();
        OrientationFull3D orientation = new OrientationFull3D(box.getSpace());
        calcOrientation(molecule, orientation);
        L.E(omega);
        rotationTensor.setOrientation(orientation);
        //find body-fixed angular momentum
        rotationTensor.transform(L);
        //now divide out moment of inertia to get body-fixed angular velocity
        omega.TE(molecule.getType().getMomentOfInertia());
        rotationTensor.invert();
        //now rotate back to get space-fixed angular velocity
        rotationTensor.transform(omega);
        return L;
    }

    @Override
    public void setAngularMomentum(IMolecule molecule, Vector com, Box box, Vector L) {
        Vector omega = ((IMoleculeOrientedKinetic)molecule).getAngularVelocity();

        RotationTensor3D rotationTensor = new RotationTensor3D();
        OrientationFull3D orientation = new OrientationFull3D(box.getSpace());
        calcOrientation(molecule, orientation);
        omega.E(L);
        rotationTensor.setOrientation(orientation);
        //find body-fixed angular momentum
        rotationTensor.transform(omega);
        //now divide out moment of inertia to get body-fixed angular velocity
        omega.DE(molecule.getType().getMomentOfInertia());
        rotationTensor.invert();
        //now rotate back to get space-fixed angular velocity
        rotationTensor.transform(omega);
    }

}
