/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.Constants;
import etomica.util.random.IRandom;

/**
 * MC move whose purpose in life is to sample an  Einstein crystal.
 * Since the energy is harmonic, each configuration can be
 * independent.
 *
 * @author Andrew Schultz
 */
public class MCMoveHO extends MCMoveBox {

    public MCMoveHO(Space space, IRandom random, double temperature, double omega) {
        super();
        this.random = random;
        this.temperature = temperature;
        this.omega2 = omega*omega;
        nBeads = box.getLeafList().size();

        lambda = new double[nBeads];
        double hbar = Constants.PLANCK_H/(2*Math.PI);
        double beta = 1.0/(Constants.BOLTZMANN_K*temperature);
        double beta_n = beta/nBeads;
        double omega_n = 1.0/(hbar*beta_n);
        double mass = box.getLeafList().get(0).getType().getMass();
        // exp(-1/2 lambda_k q_k^2)
        for(int k = 0; k < nBeads; k++){
            lambda[k] = beta_n*mass*(4.0*omega_n*omega_n*Math.sin(Math.PI*k/nBeads)*Math.sin(Math.PI*k/nBeads) + omega2);//k2=m w^2
        }

        eigenvectors = new double[nBeads][nBeads];
        for (int k = 0; k < nBeads; k++) {
            boolean doCos = k <= nBeads/2;
            for (int j = 0; j < nBeads; j++) {
                double arg = 2*Math.PI/nBeads*j*k;
                eigenvectors[k][j] = doCos ? Math.cos(arg) : -Math.sin(arg);
                eigenvectors[k][j] *= 2/Math.sqrt(nBeads);
            }
        }

    }


    public boolean doTrial() {
        double[] q = new double[nBeads];
        for (int k=0; k<nBeads; k++) {
            double fac_k = 1.0/Math.sqrt(lambda[k]);
            q[k] = fac_k*random.nextGaussian();
        }

        IAtomList atomList = box.getLeafList();
        for (int j = 0; j < nBeads; j++) {
            IAtom atomj = atomList.get(j);
            Vector rj = atomj.getPosition();
            rj.E(0);
            for (int k = 0; k < nBeads; k++) {
                rj.PE(eigenvectors[k][j]*q[k]);
            }
        }

        return true;
    }

    public double energyChange() {return 0;}

    public double getChi(double temperature) {
        return 1;
    }

    public void acceptNotify() {}

    public void rejectNotify() {}

    protected int nBeads;
    protected double temperature, omega2;
    protected final IRandom random;
    double[] lambda;
    double[][] eigenvectors;
}
