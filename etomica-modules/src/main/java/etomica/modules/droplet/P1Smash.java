/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.droplet;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.potential.IPotentialField;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * Gravity-like potential that pushes the molecules toward the center.
 */
public class P1Smash implements IPotentialField {

    public P1Smash(Space space) {
        gradient = new Vector[1];
        gradient[0] = space.makeVector();
        g = 1;
    }

    public int nBody() {
        return 1;
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

    public double virial(IAtomList atoms) {
        return 0;
    }

    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    public Vector[] gradient(IAtomList atoms) {
        IAtom a = atoms.get(0);
        if (a.getPosition().getX(2) > 0) {
            gradient[0].setX(2, g);
        }
        else {
            gradient[0].setX(2,-g);
        }
        return gradient;
    }

    public double energy(IAtomList atoms) {
        IAtom a = atoms.get(0);
        return Math.abs(a.getPosition().getX(2))*g;
    }
    
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }
    
    protected final Vector[] gradient;
    protected double g;
}
