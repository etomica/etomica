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
    public static final double hbar = 1.0; //Constants.PLANCK_H/(2.0*Math.PI);
    protected final double[] chainSigmas, dSigma;
    protected final double[] R11, R1N, dR11, dR1N;

    public MCMoveHOReal(Space space, PotentialCompute pm, IRandom random, double temperature, double omega2, Box box) {
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
        dSigma = new double[nBeads];
        R11 = new double[nBeads];
        R1N = new double[nBeads];
        dR11 = new double[nBeads];
        dR1N = new double[nBeads];

        init();
    }

    protected void init() {

        betaN = beta/nBeads;
        omegaN = 1.0/(hbar*betaN);
        double kSpring = beta * mass * omegaN*omegaN/nBeads/2;
        double dks = -kSpring/beta;
        double k0 = beta * mass * omega2/nBeads/2;
        double dk0 = k0/beta;

        double sigma2Sum = 0.5/k0;
        double dsigma2Sum = -dk0/(2*k0*k0);
        if (nBeads%2==0) {
            sigma2Sum += 0.5/(k0 + 4*kSpring);
            dsigma2Sum -= (2*dk0+8*dks)/((2*k0 + 8*kSpring)*(2*k0+8*kSpring));
        }
        for (int j=1; j*2 < nBeads; j++) {
            double s = Math.sin(Math.PI*j/nBeads);
            double d = k0 + 4*kSpring*s*s;
            sigma2Sum += 1/d;
            dsigma2Sum -= (dk0 + 4*dks*s*s)/(d*d);
        }
        sigma0 = Math.sqrt(sigma2Sum/nBeads);
        dSigma[0] = 0.5/sigma0*dsigma2Sum/nBeads;
        double bA = -Math.log(Math.sqrt(2*Math.PI)*sigma0);

        {
            chainSigmas[0] = sigma0;
            double D = -2 - k0 / kSpring;
            double dD = -dk0/kSpring + k0*dks/(kSpring*kSpring);
            double alpha = Math.log(-D/2 + Math.sqrt(D*D/4 - 1));
            double dAlpha = (-dD/2 + D*dD/(4*Math.sqrt(D*D/4-1))) / (-D/2 + Math.sqrt(D*D/4 - 1));


            // chain NM
//            System.out.println(0+" "+chainSigmas[0]+" "+dSigma[0]);
            for (int N=1; N<nBeads; N++) {
                // consider chain of length N+2 with fixed endpoints; N=n-1 corresponds to a chain from 0 back to n
                // we are computing the sigma and center coefficients for the first of the N beads
                double var1 = 0, dvar1 = 0;
                for (int l=1; l<=N; l++) {
                    double arg = l*Math.PI/(N+1);
                    double s = Math.sin(0.5*arg);
                    double lambda = k0 + 4 * kSpring * s * s;
//                    System.out.println("lambda "+k+" "+lambda[k]);
                    // vjk = sin(k*j*pi/(i+1))
                    double v0k = Math.sin(arg) / Math.sqrt(0.5 + 0.5*N);
                    double sig = v0k * Math.sqrt(0.5/lambda);
                    var1 += sig*sig;
                    dvar1 -= sig*sig / lambda * (dk0 + 4*dks*s*s);
                }
//                System.out.println(N+" "+var1+" "+dvar1);
                chainSigmas[N] = Math.sqrt(var1);
                dSigma[N] = 0.5/chainSigmas[N] * dvar1;
//                System.out.println(N+" "+chainSigmas[N]+" "+dSigma[N]);
                double d = 2*Math.sinh(alpha)*Math.sinh((N+1)*alpha);
                double dd = 2*dAlpha*Math.cosh(alpha)*Math.sinh((N+1)*alpha) + 2*(N+1)*dAlpha*Math.sinh(alpha)*Math.cosh((N+1)*alpha);
                R11[N] = (Math.cosh((N + 1)*alpha) - Math.cosh((N - 1)*alpha))/d;
                dR11[N] = ((N+1)*dAlpha*Math.sinh((N + 1)*alpha) - (N-1)*dAlpha*Math.sinh((N - 1)*alpha))/d - R11[N]/d*dd;
                R1N[N] = (Math.cosh(2*alpha) - Math.cosh(0))/d;
                dR1N[N] = (2*dAlpha*Math.sinh(2*alpha))/d - R1N[N]/d*dd;
//                System.out.println(N+" "+R1N[N]+" "+dR1N[N]);
                bA -= Math.log(Math.sqrt(2*Math.PI)*chainSigmas[N]);
            }

        }
//        System.out.println("bA: "+bA);
    }

    public double[] getChainSigmas() {
        return chainSigmas;
    }

    public double[] getDChainSigmas() {
        return dSigma;
    }

    public double[][] getCenterCoefficients() {
        return new double[][]{R11,R1N};
    }

    public double[][] getDCenterCoefficients() {
        return new double[][]{dR11,dR1N};
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
