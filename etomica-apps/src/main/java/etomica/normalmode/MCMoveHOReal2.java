/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.molecule.CenterOfMass;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculeSource;
import etomica.molecule.MoleculeSourceRandomMolecule;
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
    protected final Vector[] oldPositions;
    protected double uOld;
    protected double uaOld = Double.NaN;
    protected double uaNew = Double.NaN;
    protected double duTotal;
    protected double mass, beta, omegaN, betaN, sigma0;
    public static final double hbar = 1.0; //Constants.PLANCK_H/(2.0*Math.PI);
    protected final double[] chainSigmas, gamma, dGamma;
    protected final double[] f11, f1N, df11, df1N;
    protected final MoleculeSource moleculeSource;
    protected IMolecule molecule;
    protected Vector[] latticePositions;

    public MCMoveHOReal2(Space space, PotentialCompute pm, IRandom random, double temperature, double omega2, Box box) {
        super();
        this.pm = pm;
        this.random = random;
        this.omega2 = omega2;
        setBox(box);
        nBeads = this.box.getMoleculeList().get(0).getChildList().size();
        oldPositions = new Vector[nBeads];
        for (int j=0; j < nBeads; j++) {
            oldPositions[j] = space.makeVector();
        }

        beta = 1.0/temperature;
        // mass here is the mass of the whole ring
        mass = box.getLeafList().get(0).getType().getMass()*nBeads;

        chainSigmas = new double[nBeads];
        f11 = new double[nBeads];
        f1N = new double[nBeads];
        df11 = new double[nBeads];
        df1N = new double[nBeads];
        gamma = new double[nBeads];
        dGamma = new double[nBeads];

        moleculeSource = new MoleculeSourceRandomMolecule(box, random);

        init();

        perParticleFrequency = true;

        latticePositions = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        for (int i=0; i<latticePositions.length; i++) {
            latticePositions[i].E(CenterOfMass.position(box, box.getMoleculeList().get(i)));
        }
    }

    protected void init() {

        betaN = beta/nBeads;
        omegaN = 1.0/(hbar*betaN);
        double omegaN2 = omegaN*omegaN;

        double D = 2 + omega2 / (omegaN*omegaN);
        double alpha = Math.log(D/2 + Math.sqrt(D*D/4 - 1));
        double dAlpha = 2.0/beta/Math.sqrt(1.0+4.0*omegaN2/omega2);
        double d2Alpha = -1.0/4.0*beta*dAlpha*dAlpha*dAlpha;

        double sinhA = Math.sinh(alpha);
        double coshA = Math.cosh(alpha);
        double sinhNA = Math.sinh(nBeads*alpha);
        double coshhNA = Math.cosh(nBeads*alpha);

        double k0 = 2.0*mass*omegaN2*sinhA*Math.tanh(nBeads*alpha/2.0);
        sigma0 = k0 == 0 ? 0 : Math.sqrt(nBeads/(beta*k0));
        chainSigmas[0] = sigma0;
        gamma[0] = 1.0/2.0/beta - dAlpha/2.0*(coshA/sinhA+nBeads/sinhNA);
        dGamma[0] = -1.0/2.0/beta/beta - d2Alpha/2.0*(coshA/sinhA+nBeads/sinhNA)+ dAlpha*dAlpha/2.0*(1.0/sinhA);


        for (int i=1; i<nBeads; i++) {
            double sinhNmiA = Math.sinh((nBeads-i)*alpha);
            double coshNmiA = Math.cosh((nBeads-i)*alpha);
            double sinhNmip1A = Math.sinh((nBeads-i+1)*alpha);
            double coshNmip1A = Math.cosh((nBeads-i+1)*alpha);

            double sinhRatio = alpha == 0 ? (nBeads-i+1)/(nBeads-i) : (sinhNmip1A/sinhNmiA);
            double ki = mass*omegaN2*sinhRatio;
            chainSigmas[i] = Math.sqrt(nBeads/(beta*ki));
            gamma[i] = 1.0/2.0/beta - dAlpha/2.0*(coshNmip1A/sinhNmip1A-(nBeads-i)*sinhA/sinhNmip1A/sinhNmiA);
            f11[i] = alpha == 0 ? ((nBeads-i)/(nBeads-i+1.0)) : (sinhNmiA/sinhNmip1A);
            f1N[i] = alpha == 0 ? (1.0/(nBeads-i+1)) : (sinhA/sinhNmip1A);

            double n11 = sinhNmiA;
            double n1N = sinhA;
            double dn11 = (nBeads-i)*dAlpha*coshNmiA;
            double dn1N = dAlpha*coshA;
            double d = sinhNmip1A;
            double dd = (nBeads-i+1)*dAlpha*coshNmip1A;
            df11[i] = alpha == 0 ? 0 : ((d*dn11-n11*dd)/d/d);
            df1N[i] = alpha == 0 ? 0 : ((d*dn1N-n1N*dd)/d/d);

        }
    }

    public double[] getChainSigmas() {
        return chainSigmas;
    }

    public double[] getGamma() {
        return gamma;
    }

    public double[][] getCenterCoefficients() {
        return new double[][]{f11,f1N};
    }

    public double[][] getDCenterCoefficients() {
        return new double[][]{df11,df1N};
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
        Vector dr = box.getSpace().makeVector();
        for (int j = 0; omega2 > 0 && j < nBeads; j++) {
            Vector rj = atoms.get(j).getPosition();
            dr.Ev1Mv2(rj, latticePositions[molecule.getIndex()]);
            box.getBoundary().nearestImage(dr);
            uh += 1.0 / nBeads / 2.0 * mass * (omega2 * dr.squared());
        }
        for (int j = 0; j < nBeads; j++) {
            Vector rj = atoms.get(j).getPosition();
            int jj = j+1;
            if (jj==nBeads) jj=0;
            Vector rjj = atoms.get(jj).getPosition();
            dr.Ev1Mv2(rjj, rj);
            box.getBoundary().nearestImage(dr);
            uh += 1.0 / nBeads / 2.0 * mass * (omegaN * omegaN * dr.squared());
        }
        return uh;
    }

    public boolean doTrial() {
        molecule = moleculeSource.getMolecule();

        double uOld = pm.computeOneOldMolecule(molecule);

        IAtomList atoms = molecule.getChildList();
        Vector oldCOM = CenterOfMass.position(box, molecule);

        double uhOld = uHarmonic(molecule);
        for (int j = 0; j < nBeads; j++) {
            Vector rj = atoms.get(j).getPosition();
            oldPositions[j].E(rj);
        }

        IAtom atom0 = atoms.get(0);
        Vector prevAtomPosition = atom0.getPosition();
        for (int j = 0; j < prevAtomPosition.getD(); j++) {
            prevAtomPosition.setX(j, sigma0 * random.nextGaussian());
        }
        Vector newCOM = box.getSpace().makeVector();
        for (int k = 1; k < nBeads; k++) {
            Vector kPosition = atoms.get(k).getPosition();

            double sigma = chainSigmas[k];

            for (int j = 0; j < kPosition.getD(); j++) {
                kPosition.setX(j, sigma * random.nextGaussian());
            }

            kPosition.PEa1Tv1( f11[k], prevAtomPosition);
            kPosition.PEa1Tv1( f1N[k], atom0.getPosition());

            prevAtomPosition = kPosition;

            newCOM.PE(kPosition);
        }

        if (omega2 == 0) {
            newCOM.TE(1.0/nBeads);
            newCOM.ME(oldCOM);
            for (int k = 0; k < nBeads; k++) {
                atoms.get(k).getPosition().ME(newCOM);
            }
        }
        else {
            for (int k = 0; k < nBeads; k++) {
                atoms.get(k).getPosition().PE(latticePositions[molecule.getIndex()]);
            }
        }

        for (int k = 0; k < nBeads; k++) {
            pm.updateAtom(atoms.get(k));
        }

        double uhNew = uHarmonic(molecule);
        double uNew = pm.computeOneMolecule(molecule);
        uaOld = uOld - uhOld;
        uaNew = uNew - uhNew;
        duTotal = uNew - uOld;

        return true;
    }

    public double energyChange() {return duTotal;}

    public double getChi(double temperature) {
//        System.out.println("chi: "+Math.exp(-(uaNew - uaOld) / temperature)+" "+uaOld+" "+uaNew+" "+temperature);
        return Math.exp(-(uaNew - uaOld) / temperature);
    }

    public void acceptNotify() {
        pm.processAtomU(1);
        // put it back, then compute old contributions to energy
        IAtomList atoms = molecule.getChildList();
        Vector[] newPositions = box.getSpace().makeVectorArray(nBeads);
        for (int j=0; j<nBeads; j++) {
            Vector r = atoms.get(j).getPosition();
            newPositions[j].E(r);
            r.E(oldPositions[j]);
            pm.updateAtom(atoms.get(j));
        }
        pm.computeOneMolecule(molecule);
        pm.processAtomU(-1);
        for (int j=0; j<nBeads; j++) {
            atoms.get(j).getPosition().E(newPositions[j]);
            pm.updateAtom(atoms.get(j));
        }
    }

    public void rejectNotify() {

        IAtomList atoms = molecule.getChildList();
        for (int j=0; j<nBeads; j++) {
            atoms.get(j).getPosition().E(oldPositions[j]);
            pm.updateAtom(atoms.get(j));
        }
    }

}
