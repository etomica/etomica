/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.space.Vector;

/**
 * Atomic 3-body potential class that simply sums up contributions from multiple
 * (atomic) 3-body potentials.
 * 
 * @author Andrew Schultz
 */
public class P3AtomicSum implements Potential3Soft {

    protected final Potential3Soft[] p;

    public P3AtomicSum(Potential3Soft[] p) {
        this.p = p;
    }

    @Override
    public double u(double r212, double r213, double r223) {
        double sum = 0;
        for (int i=0; i<p.length; i++) {
            sum += p[i].u(r212, r213, r223);
        }
        return sum;
    }

    @Override
    public double u(Vector dr12, Vector dr13, Vector dr23, IAtom atom1, IAtom atom2, IAtom atom3) {
        double sum = 0;
        for (int i=0; i<p.length; i++) {
            sum += p[i].u(dr12, dr13, dr23, atom1, atom2, atom3);
        }
        return sum;
    }
}
