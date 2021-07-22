/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.atom.AtomTypeOriented;
import etomica.atom.IAtomOrientedKinetic;
import etomica.box.Box;
import etomica.space.IOrientation;
import etomica.space.Vector;

/**
 * OrientationCalc implementation that handles a monotomic oriented molecule.
 *
 * @author Andrew Schultz
 */
public class OrientationCalcAtomFasterer extends OrientationCalcNonLinear {

    private IAtomOrientedKinetic getAtom(IMolecule molecule) {
        return (IAtomOrientedKinetic)molecule.getChildList().get(0);
    }

    @Override
    public Vector getMomentOfInertia(IMolecule molecule) {
        return ((AtomTypeOriented) getAtom(molecule).getType()).getMomentOfInertia();
    }

    @Override
    public void calcOrientation(IMolecule molecule,
                                IOrientation orientation) {
        IAtomOrientedKinetic a =  getAtom(molecule);
        orientation.E(a.getOrientation());
    }

    @Override
    public void setOrientation(IMolecule molecule, Box box,
                               IOrientation orientation) {
        getAtom(molecule).getOrientation().E(orientation);
    }

    @Override
    public Vector getAngularMomentum(IMolecule molecule, Vector com, Box box) {
        Vector omega = getAtom(molecule).getAngularVelocity();
        return angularVelocityToMomentum(molecule, box, omega);
    }

    @Override
    public void setAngularMomentum(IMolecule molecule, Vector com, Box box, Vector L) {
        Vector omega = getAtom(molecule).getAngularVelocity();
        omega.E(angularMomentumToVelocity(molecule, box.getSpace(), L));
    }
}
