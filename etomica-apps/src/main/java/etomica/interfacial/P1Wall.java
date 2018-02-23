/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.interfacial;

import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.space.Space;
import etomica.space.Tensor;

public class P1Wall extends Potential1 implements PotentialSoft {

    protected final double spring, springPosition;
    protected final double gSat;
    protected final Vector[] grad;
    
    public P1Wall(Space space, double spring, double springPosition, double gSat) {
        super(space);
        this.spring = spring;
        this.springPosition = springPosition;
        this.gSat = gSat;
        grad = new Vector[]{space.makeVector()};
    }

    public double energy(IAtomList atoms) {
        double dz = atoms.get(0).getPosition().getX(2)-springPosition;
        double uSpring = 0.5*spring*dz*dz;
        return uSpring + gSat*dz;
    }

    public double virial(IAtomList atoms) {
        return 0;
    }

    public Vector[] gradient(IAtomList atoms) {
        double dz = atoms.get(0).getPosition().getX(2)-springPosition;
        grad[0].setX(2, gSat + spring*dz);
        return grad;
    }

    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

}
