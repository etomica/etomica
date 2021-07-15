/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.molecule;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.space3d.IOrientationFull3D;
import etomica.space3d.OrientationFull3D;
import etomica.space3d.RotationTensor3D;

public abstract class OrientationCalcMolecular implements OrientationCalc {

    @Override
    public Vector getAngularMomentum(IMolecule molecule, Vector com, Box box) {
        Vector dr = box.getSpace().makeVector();
        Vector L = box.getSpace().makeVector();
        // L = r x v
        for (IAtom a : molecule.getChildList()) {
            dr.Ev1Mv2(a.getPosition(), com);
            box.getBoundary().nearestImage(dr);
            dr.XE(((IAtomKinetic)a).getVelocity());
            L.PEa1Tv1(a.getType().getMass(), dr);
        }
        return L;
    }

    @Override
    public void setAngularMomentum(IMolecule molecule, Vector com, Box box, Vector L) {
        Vector omega = box.getSpace().makeVector();

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

        Vector dr = box.getSpace().makeVector();
        double omegaMag = Math.sqrt(omega.squared());
        Vector velocity = box.getSpace().makeVector();
        double mass = 0;
        for (IAtom a : molecule.getChildList()) {
            dr.Ev1Mv2(a.getPosition(), com);
            box.getBoundary().nearestImage(dr);
            double dot = dr.dot(omega);
            dr.PEa1Tv1(-dot/omegaMag, omega);
            dr.XE(omega);
            Vector v = ((IAtomKinetic)a).getVelocity();
            double m = a.getType().getMass();
            velocity.PEa1Tv1(m, v);
            mass += m;
            v.Ea1Tv1(-1, dr);
        }
        velocity.TE(1.0/mass);
        for (IAtom a : molecule.getChildList()) {
            ((IAtomKinetic)a).getVelocity().PE(velocity);
        }
    }

    @Override
    public void setOrientation(IMolecule molecule, Box box,
                               IOrientationFull3D orientation) {
        Vector com = CenterOfMass.position(box, molecule);
        IAtomList childList = molecule.getChildList();
        molecule.getType().getConformation().initializePositions(childList);
        Vector com0 = CenterOfMass.position(box, molecule);
        RotationTensor3D rotationTensor = new RotationTensor3D();
        rotationTensor.setOrientation(orientation);
        rotationTensor.invert();
        for (int iChild = 0; iChild<childList.size(); iChild++) {
            IAtom a = childList.get(iChild);
            Vector r = a.getPosition();
            r.ME(com0);
            rotationTensor.transform(r);
            r.PE(com);
        }
    }
}
