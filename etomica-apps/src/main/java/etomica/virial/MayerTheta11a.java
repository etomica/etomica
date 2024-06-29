/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotentialMolecular;
import etomica.potential.compute.PotentialCompute;

/**
 * MayerFunction implementation that returns f * du/dk where
 * u is intramolecular energy and k is some intramolecular parameter.
 * Optionally, U*exp(-beta*U) can be added, where U is the intermolecular energy.
 *
 * This class is useful when computing dbeta/dk at fixed B2 (B2=0 for theta point).
 */
public class MayerTheta11a implements MayerFunction {

    protected final IPotentialMolecular potential;
    // pcdk returns du/dk for the intramolecular energy of all molecules, probably a PotentialMasterBonding
    // pcu returns intramolecular energy
    protected PotentialCompute pcdk, pcu;

    /**
     * Constructor Mayer function using given potential.
     *   addUE is true to get dB2/dbeta
     *   addUE is false to get dB2/dk
     */
    public MayerTheta11a(IPotentialMolecular potential) {
        this.potential = potential;
    }

    public void setPotentialDK(PotentialCompute pcdk) {
        this.pcdk = pcdk;
    }
    public void setPotentialu(PotentialCompute pcu) {
        this.pcu = pcu;
    }

    // e*beta*U*dudk + f*dudk*(beta*u-1)
    public double f(IMoleculeList pair, double r2, double beta) {
        double U = potential.energy(pair);
        double x = -beta*U;
        double f;
        if (Math.abs(x) < 0.01) {
            f = x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;
        }
        else {
            f = Math.exp(x) - 1;
        }
        if ((Double.isNaN(f) || Double.isInfinite(f))) {
            throw new  RuntimeException("bogus f: "+f+"   beta: "+beta+"   u: "+potential.energy(pair));
        }
        double e = f+1;
        double dudk = pcdk.computeAll(false);
        double u = pcu.computeAll(false);
        return -e*x*dudk + f*dudk*(beta*u-1);
    }

    public void setBox(Box newBox) {
    }
}
