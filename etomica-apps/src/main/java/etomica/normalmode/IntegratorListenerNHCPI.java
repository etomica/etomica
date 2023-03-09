/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.box.Box;
import etomica.integrator.IntegratorListenerNHC;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class IntegratorListenerNHCPI extends IntegratorListenerNHC {

    /**
     * Creates a Nose-Hoover chains thermostat
     *
     * @param integrator
     * @param random
     * @param numChains  # of chains in the thermostat
     * @param tau        thermostat timescale
     */
    public IntegratorListenerNHCPI(IntegratorPIMD integrator, IRandom random, int numChains, double tau) {
        super(integrator, random, numChains, tau);
    }

    protected void computeKE() {
        Box box = integrator.getBox();
        double[] mScale = ((IntegratorPIMD)integrator).getMassScale();
        //printf("fac %e  %e %e\n", fac, etaP[0], q[0]);
        kineticEnergy = 0;
        int n = 0;
        int D = box.getSpace().getD();

        for (IAtom atom : box.getLeafList()) {
            IAtomKinetic a = (IAtomKinetic) atom;
            double mass = a.getType().getMass() * mScale[a.getIndex()];
            if (mass == Double.POSITIVE_INFINITY) continue;
            // rigid molecules will confuse us here
            n += D;
            Vector v = a.getVelocity();
            kineticEnergy += 0.5 * mass * v.squared();
        }
        if (n != ndof) {
            ndof = n;
            double x = integrator.getTemperature() * tau * tau;
            q[0] = ndof * x;
        }
    }

    protected void propagatorU3(double dt) {
        Box box = integrator.getBox();
        double[] mScale = ((IntegratorPIMD)integrator).getMassScale();
        double fac = Math.exp(-dt * etaP[0] / q[0]);
        //printf("fac %e  %e %e\n", fac, etaP[0], q[0]);
        kineticEnergy = 0;
        for (IAtom atom : box.getLeafList()) {
            IAtomKinetic a = (IAtomKinetic) atom;
            double mass = a.getType().getMass() * mScale[a.getIndex()];
            if (mass == Double.POSITIVE_INFINITY) continue;
            Vector v = a.getVelocity();
            v.TE(fac);
            kineticEnergy += 0.5 * mass * v.squared();
        }
        for (int i = 0; i < q.length; i++) {
            eta[i] += dt * etaP[i] / q[i];
        }
    }

}
