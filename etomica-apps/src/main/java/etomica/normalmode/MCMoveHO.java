/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class MCMoveHO extends MCMoveBox {
    protected int nBeads;
    protected double temperature, omega2;
    protected final IRandom random;
    double[] lambdaN;
    double[][] eigenvectors;
    PotentialCompute pm;
    protected final Vector[] oldPositions;
    protected double uOld;
    protected double uaOld = Double.NaN;
    protected double uaNew = Double.NaN;
    protected double mass, beta, omegaN, betaN;
    protected Box box;
    public static final double hbar = 1.0; //Constants.PLANCK_H/(2.0*Math.PI);
    public static final double kB = 1.0;//Constants.BOLTZMANN_K


    public MCMoveHO(Space space, PotentialCompute pm, IRandom random, double temperature, double omega, Box box) {
        super();
        this.pm = pm;
        this.random = random;
        this.temperature = temperature;
        this.omega2 = omega*omega;
        this.box = box;
        nBeads = this.box.getMoleculeList().get(0).getChildList().size();
        oldPositions = new Vector[nBeads];
        for (int i=0; i < nBeads; i++){
            oldPositions[i] = space.makeVector();
        }

        lambdaN = new double[nBeads];
        beta = 1.0/(kB*temperature);
        betaN = beta/nBeads;
        omegaN = 1.0/(hbar*betaN);
        mass = box.getLeafList().get(0).getType().getMass();
        // exp(- betan * 1/2*lambdaN_k q_k^2)
        int nK = nBeads/2;
        lambdaN[nK] = mass*omega2; //k=0
        for(int k = 1; k <= nK; k++){ //-2...2
            double lambda_k = mass*(4.0*omegaN*omegaN*Math.sin(Math.PI*k/nBeads)*Math.sin(Math.PI*k/nBeads) + omega2);
            lambdaN[nK-k] = lambda_k;
            if (k != nK || nBeads % 2 != 0){ //odd
                lambdaN[nK+k] = lambda_k;
            }

        }

        eigenvectors = new double[nBeads][nBeads];
        for (int i = 0; i < nBeads; i++) {
            eigenvectors[i][nK] = 1.0/Math.sqrt(nBeads);//k=0
            for (int k = 1; k <= nK; k++) {
                double arg = 2.0*Math.PI/nBeads*i*k;
                eigenvectors[i][nK-k] =  2*Math.sin(-arg)/Math.sqrt(nBeads);
                if (k != nK || nBeads % 2 != 0){ //odd
                    eigenvectors[i][nK+k] =  2*Math.cos(arg)/Math.sqrt(nBeads);
                }
            }
        }
    }


    public boolean doTrial() {
        double uOld = pm.getLastEnergy();
        IAtomList atomList = box.getLeafList();
        IAtom atomi, atomip1;
        Vector ri, rip1;
        double uhOld = 0;
        for (int i = 0; i < nBeads; i++) {
            atomi = atomList.get(i);
            atomip1 = (i == nBeads - 1) ? atomList.get(0) : atomList.get(i + 1);
            ri = atomi.getPosition();
            oldPositions[i].E(ri);
            rip1 = atomip1.getPosition();
            Vector dr = box.getSpace().makeVector();
            dr.Ev1Mv2(ri, rip1);
            uhOld += 1.0 / nBeads / 2.0 * mass * (omegaN * omegaN * dr.squared() + omega2 * ri.squared());
        }
        uaOld = uOld - uhOld;

        double uhNew = 0;
        double[] q = new double[nBeads];
        int nK = nBeads/2;
        double fack = 1.0/Math.sqrt(betaN * lambdaN[nK]);
        q[nK] = fack*random.nextGaussian();
        uhNew += 1.0/2.0/beta*betaN*lambdaN[nK]*(q[nK]*q[nK]);
        for (int k=1; k<=nK; k++) {
            fack = 1.0/Math.sqrt(2*betaN * lambdaN[nK-k]);
            q[nK-k] = fack*random.nextGaussian();
            uhNew += 1.0/beta*betaN*lambdaN[nK-k]*q[nK-k]*q[nK-k];
            if (k != nK || nBeads % 2 != 0){ //odd
                fack = 1.0/Math.sqrt(2*betaN * lambdaN[nK+k]);
                q[nK+k] = fack*random.nextGaussian();
                uhNew += 1.0/beta*betaN*lambdaN[nK+k]*q[nK+k]*q[nK+k];
            }
        }

        for (int i = 0; i < nBeads; i++) {
            atomi = atomList.get(i);
            ri = atomi.getPosition();
            ri.E(0);
            for (int k = 0; k < nBeads; k++) {
                ri.PE(eigenvectors[i][k]*q[k]);
            }
        }
        pm.init();
        double uNew = pm.computeAll(false);
        uaNew = uNew - uhNew;

        return true;
    }

    public double energyChange() {return 0;}

    public double getChi(double temperature) {
        return Math.exp(-(uaNew - uaOld) / (kB*temperature));

    }

    public void acceptNotify() { /* do nothing */}

    public void rejectNotify() {
        IAtomList atomList = box.getLeafList();
        IAtom atomi;
        for (int i = 0; i < nBeads; i++) {
            atomi = atomList.get(i);
            atomi.getPosition().E(oldPositions[i]);
        }
        pm.init();
        pm.computeAll(false);
    }

}
