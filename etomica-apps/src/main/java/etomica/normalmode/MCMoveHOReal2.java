/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Move that generates coordinates for PI rings in an external field.  Coordinates of the atoms are generated
 * sequentially, each from a Gaussian distribution.  The cost is O(n), where n is the number of beads in the ring.
 * The move is accepted/rejected based on the difference between the actual energy of the system and the harmonic
 * approximation to the energy used to generate the configuration.
 */
public class MCMoveHOReal2 extends MCMoveBox {
    protected int nBeads;
    protected double omega2;
    protected final IRandom random;
    PotentialCompute pm;
    protected final Vector[][] oldPositions;
    protected double uOld;
    protected double uaOld = Double.NaN;
    protected double uaNew = Double.NaN;
    protected double duTotal;
    protected double mass, beta, omegaN, betaN, sigma0;
    public static final double hbar = 1.0; //Constants.PLANCK_H/(2.0*Math.PI);
    protected final double[] chainSigmas;
    protected final double[] f11, f1N;

    public MCMoveHOReal2(Space space, PotentialCompute pm, IRandom random, double temperature, double omega2, Box box) {
        super();
        this.pm = pm;
        this.random = random;
        this.omega2 = omega2;
        setBox(box);
        nBeads = this.box.getMoleculeList().get(0).getChildList().size();
        oldPositions = new Vector[box.getMoleculeList().size()][nBeads];
        for (int i=0; i<oldPositions.length; i++) {
            for (int j=0; j < nBeads; j++) {
                oldPositions[i][j] = space.makeVector();
            }
        }

        beta = 1.0/temperature;
        mass = box.getLeafList().get(0).getType().getMass();

        chainSigmas = new double[nBeads];
        f11 = new double[nBeads];
        f1N = new double[nBeads];

        init();
    }

    protected void init() {

        betaN = beta/nBeads;
        omegaN = 1.0/(hbar*betaN);

        double D = 2 + omega2 / (omegaN*omegaN);
        double alpha = Math.log(D/2 + Math.sqrt(D*D/4 - 1));
        double C0 = mass*(omega2 + 2*omegaN*omegaN)*Math.tanh(alpha)*Math.tanh(nBeads*alpha/2);
        sigma0 = C0 == 0 ? 0 : Math.sqrt(nBeads/(beta*C0));
        chainSigmas[0] = sigma0;

        for (int i=1; i<nBeads; i++) {
            double tanhRatio = alpha == 0 ? 1.0/(nBeads-i) : (Math.tanh(alpha)/Math.tanh((nBeads-i)*alpha));
            double Ci = 0.5*mass*(omega2 + 2*omegaN*omegaN)*(1 + tanhRatio);
            chainSigmas[i] = Math.sqrt(nBeads/(beta*Ci));

            f11[i] = alpha == 0 ? ((nBeads-i)/(nBeads-i+1.0)) : (Math.sinh((nBeads-i)*alpha)/Math.sinh((nBeads-i+1)*alpha));
            f1N[i] = alpha == 0 ? (1.0/(nBeads-i+1)) : (Math.sinh(alpha)/Math.sinh((nBeads-i+1)*alpha));
        }
    }

    public double[] getChainSigmas() {
        return chainSigmas;
    }

    public double[][] getCenterCoefficients() {
        return new double[][]{f11,f1N};
    }

    public void setOmega2(double omega2) {
        this.omega2 = omega2;
        init();
    }

    public void setTemperature(double temperature) {
        this.beta = 1/temperature;
        init();
    }

    protected double uHarmonic(IMolecule molecule) {
        IAtomList atoms = molecule.getChildList();
        double uh = 0;
        for (int j = 0; j < nBeads; j++) {
            Vector rj = atoms.get(j).getPosition();
            uh += 1.0 / nBeads / 2.0 * mass * (omega2 * rj.squared());
        }
        for (int j = 0; j < nBeads; j++) {
            Vector rj = atoms.get(j).getPosition();
            int jj = j+1;
            if (jj==nBeads) jj=0;
            Vector rjj = atoms.get(jj).getPosition();
            uh += 1.0 / nBeads / 2.0 * mass * (omegaN * omegaN * rjj.Mv1Squared(rj));
        }
        return uh;
    }

    public boolean doTrial() {
        double uOld = pm.getLastEnergy();
        IMoleculeList molecules = box.getMoleculeList();

        double uhOld = 0, uhNew = 0;

        for (int i = 0; i<molecules.size(); i++) {
            IAtomList atoms = molecules.get(i).getChildList();
            Vector COM = box.getSpace().makeVector();

            uhOld += uHarmonic(molecules.get(i));
            for (int j = 0; j < nBeads; j++) {
                Vector rj = atoms.get(j).getPosition();
                oldPositions[i][j].E(rj);
                COM.PE(rj);
            }

            IAtom atom0 = atoms.get(0);
            oldPositions[i][0].E(atom0.getPosition());
            Vector prevAtomPosition = atom0.getPosition();
            for (int j = 0; j < prevAtomPosition.getD(); j++) {
                prevAtomPosition.setX(j, sigma0 * random.nextGaussian());
            }
            for (int k = 1; k < nBeads; k++) {
                Vector kPosition = atoms.get(k).getPosition();
                oldPositions[i][k].E(kPosition);

                double sigma = chainSigmas[k];

                for (int j = 0; j < kPosition.getD(); j++) {
                    kPosition.setX(j, sigma * random.nextGaussian());
                }

                kPosition.PEa1Tv1( f11[k], prevAtomPosition);
                kPosition.PEa1Tv1( f1N[k], atom0.getPosition());

                prevAtomPosition = kPosition;

                COM.ME(kPosition);
            }
            uhNew += uHarmonic(molecules.get(i));
            COM.TE(1.0/nBeads);

            if (omega2 == 0) {
                for (int k = 0; k < nBeads; k++) {
                    atoms.get(k).getPosition().PE(COM);
                }
            }
        }
        pm.init();
        double uNew = pm.computeAll(false);
        uaOld = uOld - uhOld;
        uaNew = uNew - uhNew;
        duTotal = uNew - uOld;

        return true;
    }

    public double energyChange() {return duTotal;}

    public double getChi(double temperature) {
        return Math.exp(-(uaNew - uaOld) / temperature);

    }

    public void acceptNotify() { /* do nothing */}

    public void rejectNotify() {
        IMoleculeList molecules = box.getMoleculeList();
        for (int i = 0; i<molecules.size(); i++) {
            IAtomList atoms = molecules.get(i).getChildList();
            for (int j=0; j<nBeads; j++) {
                atoms.get(j).getPosition().E(oldPositions[i][j]);
            }
        }
        pm.init();
        pm.computeAll(false);
    }

}
