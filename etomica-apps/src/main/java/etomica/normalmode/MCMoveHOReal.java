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
public class MCMoveHOReal extends MCMoveBox {
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
    protected Box box;
    public static final double hbar = 1.0; //Constants.PLANCK_H/(2.0*Math.PI);
    protected final double[] chainSigmas;
    protected final double[] R11, R1N;

    public MCMoveHOReal(Space space, PotentialCompute pm, IRandom random, double temperature, double omega2, Box box) {
        super();
        this.pm = pm;
        this.random = random;
        this.omega2 = omega2;
        this.box = box;
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
        R11 = new double[nBeads];
        R1N = new double[nBeads];

        init();
    }

    protected void init() {

        betaN = beta/nBeads;
        omegaN = 1.0/(hbar*betaN);

        {
            double[] lambdaN = new double[nBeads];
            // exp(- betan * 1/2*lambdaN_k q_k^2)
            for(int k = 0; k <= nBeads/2; k++){ //-2...2
                lambdaN[k] = mass*(4.0*omegaN*omegaN*Math.sin(Math.PI*k/nBeads)*Math.sin(Math.PI*k/nBeads) + omega2);
//                System.out.println("real move "+k+" "+lambdaN[k]);
            }

            // eigenvector component for atom 0 is always 1/sqrt(n)
            double evec0 = 1 / Math.sqrt(nBeads);

            // for a given atom i and wave vector, sigma_i = ev_i * sigma
            double sigma2Sum = 0;
            for (int k = 0; k <= nBeads/2; k++) {
                double onetwo = k==0 || 2*k==nBeads ? 1 : 2;
                double fack = Math.sqrt(onetwo/(betaN*lambdaN[k]));
                double sigmai = fack * evec0;
                sigma2Sum += sigmai*sigmai;
            }
            sigma0 = Math.sqrt(sigma2Sum);
//            System.out.println("sigma0 "+sigma0+" keff "+0.5/(sigma0*sigma0));
        }

        {
            chainSigmas[0] = sigma0;
            double kSpring = beta * mass * omegaN*omegaN/nBeads/2;
            double k0 = beta * mass * omega2/nBeads/2;
            double D = -2 - k0 / kSpring;
            double alpha = Math.log(-D/2 + Math.sqrt(D*D/4 - 1));

            // chain NM
            for (int i=1; i<nBeads; i++) {
                // consider chain of length i+2 with fixed endpoints; i=n-1 corresponds to a chain from 0 to n
                double var1 = 0;
                for (int k=0; k<i; k++) {
                    double s = Math.sin((k + 1) * Math.PI / (2 * (i + 1)));
                    double lambdak = k0 + 4 * kSpring * s * s;
//                    System.out.println("lambda "+k+" "+lambda[k]);
                    // vjk = sin(k*j*pi/(i+1))
                    double v0k = Math.sin((k+1)*(0+1)*Math.PI/(i+1)) / Math.sqrt(0.5 + 0.5*i);
                    double sig = v0k * Math.sqrt(0.5/lambdak);
                    var1 += sig*sig;
                }
                chainSigmas[i] = Math.sqrt(var1);
//                System.out.println("chain sigma "+i+" "+chainSigmas[i]);
                double d = 2*Math.sinh(alpha)*Math.sinh((i+1)*alpha);
                R11[i] = (Math.cosh((i + 1)*alpha) - Math.cosh((i - 1)*alpha))/d;
                R1N[i] = (Math.cosh(2*alpha) - Math.cosh(0))/d;
            }
//            System.out.println("D "+D+"   alpha "+alpha);

        }
    }

    public double[] getChainSigmas() {
        return chainSigmas;
    }

    public double[][] getCenterCoefficients() {
        return new double[][]{R11,R1N};
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

            uhOld += uHarmonic(molecules.get(i));
            for (int j = 0; j < nBeads; j++) {
                Vector rj = atoms.get(j).getPosition();
                oldPositions[i][j].E(rj);
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

                int N = nBeads - k; // # of beads that still need to be inserted
                double sigma = chainSigmas[N];

                for (int j = 0; j < kPosition.getD(); j++) {
                    kPosition.setX(j, sigma * random.nextGaussian());
                }

                kPosition.PEa1Tv1( R11[N], prevAtomPosition);
                kPosition.PEa1Tv1( R1N[N], atom0.getPosition());

                prevAtomPosition = kPosition;
            }
            uhNew += uHarmonic(molecules.get(i));
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
