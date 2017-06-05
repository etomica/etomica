/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * Hard potential class that wraps another hard potential.
 *
 * @author Andrew Schultz
 */
 public class P2HardWrapper implements PotentialHard {

    public P2HardWrapper(Space space, PotentialHard potential) {
        this.space = space;
        wrappedPotential = potential;
    }

    public double energy(IAtomList atoms) {
        return wrappedPotential.energy(atoms);
    }

    public int nBody() {
        return wrappedPotential.nBody();
    }

    public void setBox(Box box) {
        wrappedPotential.setBox(box);
    }
    
    public PotentialHard getWrappedPotential() {
        return wrappedPotential;
    }

    public void setWrappedPotential(PotentialHard newWrappedPotential) {
        wrappedPotential = newWrappedPotential;
    }

    public double getRange() {
        return wrappedPotential.getRange();
    }
    
    public void bump(IAtomList atoms, double falseTime) {
        wrappedPotential.bump(atoms, falseTime);
    }

    public double collisionTime(IAtomList atoms, double falseTime) {
        return wrappedPotential.collisionTime(atoms, falseTime);
    }

    public double energyChange() {
        return wrappedPotential.energyChange();
    }

    public double lastCollisionVirial() {
        return wrappedPotential.lastCollisionVirial();
    }

    public Tensor lastCollisionVirialTensor() {
        return wrappedPotential.lastCollisionVirialTensor();
    }

    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected PotentialHard wrappedPotential;
}
