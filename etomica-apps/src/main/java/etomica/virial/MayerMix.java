/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotentialMolecular;

/**
 * @author kofke
 *
 * Mayer function that invokes different potentials for different species
 * molecule pairs (based on the species).
 */
public class MayerMix implements MayerFunction {

    private final IPotentialMolecular[][] potential;
    protected static final boolean debug = false;

    /**
     * Constructor Mayer function using given potential.
     */
    public MayerMix(IPotentialMolecular[][] potential) {
        this.potential = potential;
    }

    public double f(IMoleculeList pair, double r2, double beta) {
        int idx0 = pair.get(0).getType().getIndex();
        int idx1 = pair.get(1).getType().getIndex();
        double x = -beta*potential[idx0][idx1].energy(pair);
        double f;
        if (Math.abs(x) < 0.01) {
            f = x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;
        }
        else {
            f = Math.exp(x) - 1;
        }
        if (debug && (Double.isNaN(f) || Double.isInfinite(f))) {
            throw new  RuntimeException("bogus f: "+f+"   beta: "+beta+"   u: "+potential[idx0][idx1].energy(pair));
        }
        return f;
    }

    public void setBox(Box newBox) {
    }
}
