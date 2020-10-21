/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.atom.AtomTest;
import etomica.atom.IAtom;
import etomica.util.random.IRandom;

/**
 * This class accepts any atom with a random probability.
 */
public class AtomTestRandom implements AtomTest {

    protected final IRandom random;
    protected double frac;

    public AtomTestRandom(IRandom random) {
        this.random = random;
    }

    public void setFraction(double newFrac) {
        this.frac = newFrac;
    }

    public double getFraction() {
        return frac;
    }

    @Override
    public boolean test(IAtom a) {
        return frac > random.nextDouble();
    }
}
