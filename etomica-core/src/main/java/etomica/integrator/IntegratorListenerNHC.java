/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.integrator;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.space.Vector;
import etomica.units.dimensions.Energy;
import etomica.util.random.IRandom;

/**
 * Listener that handles Nose-Hoover chains thermostat.
 * adapted from Allen & Tildesley's example code, md_nvt_lj
 * https://github.com/Allen-Tildesley/examples
 *
 * based on Martyna et al, Molec Phys, 87, 1117 (1996)  https://doi.org/10.1080/00268979600100761
 * and Tuckerman et al J Phys A, 39, 5629 (2006)        https://doi.org/10.1088/0305-4470/39/19/S18
 */
public class IntegratorListenerNHC implements IntegratorListenerMD {

    protected final IntegratorMD integrator;
    protected final IRandom random;
    protected double kineticEnergy;
    protected double[] q;         // thermal inertias for each chain
    protected double[] eta;       // thermal coordinates
    protected double[] etaP;      // thermal momenta
    protected double tau;         // thermostat timescale
    protected int ndof;
    protected double chainEnergy;

    /**
     * Creates a Nose-Hoover chains thermostat
     *
     * @param numChains  # of chains in the thermostat
     * @param tau        thermostat timescale
     */
    public IntegratorListenerNHC(IntegratorMD integrator, IRandom random, int numChains, double tau) {
        this.integrator = integrator;
        this.random = random;
        this.tau = tau;
        setNumChains(numChains);
    }

    public void setTau(double tau) {
        this.tau = tau;
        reset();
    }

    public void setNumChains(int numChains) {
        q = new double[numChains];
        eta = new double[numChains];
        etaP = new double[numChains];
        reset();
    }

    protected void reset() {
        double temperature = integrator.getTemperature();
        double x = temperature * tau * tau;
        q[0] = ndof * x;
        for (int i = 1; i < q.length; i++) q[i] = x;
        double sqrtT = Math.sqrt(temperature);
        for (int i = 0; i < etaP.length; i++) {
            eta[i] = 0;
            etaP[i] = random.nextGaussian() * sqrtT;
        }
    }

    public double[] getEtaP() {
        return etaP;
    }

    public double[] getEta() {
        return eta;
    }

    public void setEtaP(double[] etaP) {
        System.arraycopy(etaP, 0, this.etaP, 0, etaP.length);
    }

    public double getChainEnergy() {
        return chainEnergy;
    }

    @Override
    public void integratorForcePrecomputed(IntegratorEvent e) {
        doStuff();
    }

    @Override
    public void preThermostat(IntegratorEvent e) {
        doStuff();

        // PE + KE + CE = constant (should be conserved)
        chainEnergy = 0;
        for (int i=0; i<q.length; i++) {
            chainEnergy += 0.5*etaP[i]*etaP[i]/q[i];
        }
        double sum2 = ndof * eta[0];
        for (int i=1; i<eta.length; i++) {
            sum2 += eta[i];
        }
        chainEnergy += integrator.getTemperature() * sum2;
    }

    protected void doStuff() {
        double dt = integrator.getTimeStep();
        // superclass integrator only computes KE after step
        computeKE();
        propagatorU4(dt / 4, -1);
        propagatorU3(dt / 2);
        propagatorU4(dt / 4, +1);
    }

    protected void computeKE() {
        Box box = integrator.getBox();
        //printf("fac %e  %e %e\n", fac, etaP[0], q[0]);
        kineticEnergy = 0;
        int n = 0;
        int D = box.getSpace().getD();
        for (IAtom atom : box.getLeafList()) {
            IAtomKinetic a = (IAtomKinetic) atom;
            double mass = a.getType().getMass();
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
        double fac = Math.exp(-dt * etaP[0] / q[0]);
        //printf("fac %e  %e %e\n", fac, etaP[0], q[0]);
        kineticEnergy = 0;
        for (IAtom atom : box.getLeafList()) {
            IAtomKinetic a = (IAtomKinetic) atom;
            double mass = a.getType().getMass();
            if (mass == Double.POSITIVE_INFINITY) continue;
            Vector v = a.getVelocity();
            v.TE(fac);
            kineticEnergy += 0.5 * mass * v.squared();
        }
        for (int i = 0; i < q.length; i++) {
            eta[i] += dt * etaP[i] / q[i];
        }
    }

    protected void propagatorU4(double dt, int direction) {
        int numChains = q.length;
        double temperature = integrator.getTemperature();
        for (int i = 0; i < numChains; i++) {
            int ic = direction == 1 ? i : (numChains - i - 1);
            double gi = ic == 0 ? (2 * kineticEnergy - ndof * temperature) : ((etaP[ic - 1] * etaP[ic - 1] / q[ic - 1]) - temperature);
            if (ic == numChains - 1) {
                etaP[ic] += dt * gi;
                //printf("etaP[0] %e\n", etaP[0]);
                continue;
            }
            double x = dt * etaP[ic + 1] / q[ic + 1];
            double c;
            double expmx = Math.exp(-x);
            if (x < 0.001) {
                // for small values, use Taylor series
                c = 1 + x * (-0.5 + x * (1.0 / 6.0 - x / 24));
            } else {
                c = (1 - expmx) / x;
            }
            etaP[ic] = etaP[ic] * expmx + dt * gi * c;
        }
    }

    /**
     * DataSource for the NHC external energy.
     */
    public static class DataSourceChainEnergy extends DataSourceScalar {

        protected final IntegratorListenerNHC nhc;

        public DataSourceChainEnergy(IntegratorListenerNHC nhc) {
            super("NHC Energy", Energy.DIMENSION);
            this.nhc = nhc;
        }

        @Override
        public double getDataAsScalar() {
            return nhc.getChainEnergy();
        }
    }

    public static class DataSourceTotalEnergy extends DataSourceScalar {

        protected final IntegratorMD integrator;
        protected final IntegratorListenerNHC nhc;

        public DataSourceTotalEnergy(IntegratorMD integrator, IntegratorListenerNHC nhc) {
            super("Total Energy", Energy.DIMENSION);
            this.integrator = integrator;
            this.nhc = nhc;
        }

        @Override
        public double getDataAsScalar() {
            return integrator.getPotentialEnergy() + integrator.getKineticEnergy() + nhc.getChainEnergy();
        }
    }
}
