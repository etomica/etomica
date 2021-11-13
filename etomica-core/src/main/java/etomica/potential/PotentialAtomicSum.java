/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;

/**
 * Atomic potential class that simply sums up contributions from multiple
 * (atomic) potentials.
 * 
 * @author Andrew Schultz
 */
public class PotentialAtomicSum implements IPotentialAtomic {

    protected final IPotentialAtomic[] p;
    
    public PotentialAtomicSum(IPotentialAtomic[] p) {
        this.p = p;
    }

    public double getRange() {
        double r = 0;
        for (int i=0; i<p.length; i++) {
            if (r < p[i].getRange()) r = p[i].getRange();
        }
        return r;
    }

    public int nBody() {
        return p[0].nBody();
    }

    public double energy(IAtomList atoms) {
        double sum = 0;
        for (int i=0; i<p.length; i++) {
            sum += p[i].energy(atoms);
        }
        return sum;
    }

}
