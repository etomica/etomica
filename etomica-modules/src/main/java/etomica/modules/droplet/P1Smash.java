/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.droplet;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.potential.IPotentialField;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Gravity-like potential that pushes the molecules toward the center.
 */
public class P1Smash implements IPotentialField {

    public P1Smash(Space space) {
        g = 1;
    }

    public void setG(double newG) {
        g = newG;
    }

    public double getG() {
        return g;
    }

    public double u(IAtom atom) {
        return Math.abs(atom.getPosition().getX(2)) * g;
    }

    /**
     * Computes the force (and adds it to f) for IAtom atom and returns the
     * energy due to the field.
     */
    public double udu(IAtom atom, Vector f) {
        double z = atom.getPosition().getX(2);
        f.setX(2, f.getX(2) - g * Math.signum(z));
        return Math.abs(z) * g;
    }

    public double energy(IAtomList atoms) {
        IAtom a = atoms.get(0);
        return Math.abs(a.getPosition().getX(2))*g;
    }

    protected double g;
}
