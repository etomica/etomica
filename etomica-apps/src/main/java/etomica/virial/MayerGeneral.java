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
 * General Mayer function, which wraps the Mayer potential around an instance of
 * a Potential2 object.
 */
public class MayerGeneral implements MayerFunction {

    protected static final boolean debug = false;
    
    /**
     * Constructor Mayer function using given potential.
     */
    public MayerGeneral(IPotentialMolecular potential) {
        this.potential = potential;
    }

    public double f(IMoleculeList pair, double r2, double beta) {
        double x = -beta*potential.energy(pair);
        double f;
        if (Math.abs(x) < 0.01) {
            f = x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;
        }
        else {
            f = Math.exp(x) - 1;
        }
        if (debug && (Double.isNaN(f) || Double.isInfinite(f))) {
            throw new  RuntimeException("bogus f: "+f+"   beta: "+beta+"   u: "+potential.energy(pair));
        }
        return f;
    }

    public void setBox(Box newBox) {
    }

    private final IPotentialMolecular potential;
}
