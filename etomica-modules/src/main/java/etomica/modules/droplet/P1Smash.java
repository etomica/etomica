/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.droplet;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.potential.PotentialSoft;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * Gravity-like potential that pushes the molecules toward the center.
 */
public class P1Smash implements PotentialSoft {

    public P1Smash(Space space) {
        gradient = new Vector[1];
        gradient[0] = space.makeVector();
        g = 1;
    }
    
    public void setBox(Box newBox) {}
    
    public int nBody() {
        return 1;
    }
    
    public void setG(double newG) {
        g = newG;
    }
    
    public double getG() {
        return g;
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
