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
    protected double hbar;
    protected double mass, beta, omegaN, sigma0;
    protected final double[] chainSigmas, gamma, dGamma;
    protected final double[] f11, f1N, df11, df1N,d2f11, d2f1N;
    protected final MoleculeSource moleculeSource;
    protected IMolecule molecule;
    protected Vector[] latticePositions;
    protected int nGrow;

    public MCMoveHOReal2(Space space, PotentialCompute pm, IRandom random, double temperature, double omega2, Box box, double hbar) {
        super();
        this.pm = pm;
        this.random = random;
        this.omega2 = omega2;
        this.hbar = hbar;
        setBox(box);
        nBeads = this.box.getMoleculeList().get(0).getChildList().size();
        oldPositions = new Vector[nBeads];
        for (int j=0; j < nBeads; j++) {
            oldPositions[j] = space.makeVector();
        }

        beta = 1.0/temperature;
        // mass here is the mass of the whole ring
        mass = nBeads * box.getLeafList().get(0).getType().getMass();  // ring mass, m

        chainSigmas = new double[nBeads];
        f11 = new double[nBeads];
        f1N = new double[nBeads];
        df11 = new double[nBeads];
        df1N = new double[nBeads];
        d2f11 = new double[nBeads];
        d2f1N = new double[nBeads];
        gamma = new double[nBeads];
        dGamma = new double[nBeads];

        moleculeSource = new MoleculeSourceRandomMolecule(box, random);

        init();

        perParticleFrequency = true;

        latticePositions = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        for (int i=0; i<latticePositions.length; i++) {
            latticePositions[i].E(CenterOfMass.position(box, box.getMoleculeList().get(i)));
        }

        nGrow = 1;
        if (nBeads > 1) {
            nGrow = nBeads;
        }
        System.out.println(" nGrow: " + nGrow);
    }

    public void setNumGrow(int nGrow) {
        this.nGrow = nGrow;
    }

    protected void init() {
        omegaN = Math.sqrt(nBeads)/(hbar*beta);
        double D = 2 + omega2 / (nBeads*omegaN*omegaN);
        double alpha = Math.log(D/2 + Math.sqrt(D*D/4 - 1));
        double dAlpha = 2.0/beta/Math.sqrt(1.0+4.0*nBeads*omegaN*omegaN/omega2);
        double dAlpha2 = dAlpha*dAlpha;
        double d2Alpha = -1.0/4.0*beta*dAlpha*dAlpha*dAlpha;

        double sinhA = Math.sinh(alpha);
        double coshA = Math.cosh(alpha);
        double sinhNA = Math.sinh(nBeads*alpha);
        double coshhNA = Math.cosh(nBeads*alpha);

        double k0 = 2*mass*omegaN*omegaN*sinhA*Math.tanh(nBeads*alpha/2.0);
        sigma0 = k0 == 0 ? 0 : 1/Math.sqrt(beta*k0);
        chainSigmas[0] = sigma0;

        gamma[0] = alpha == 0 ? 0 : 1.0/2.0/beta - dAlpha/2.0*(coshA/sinhA+nBeads/sinhNA);
        dGamma[0] = alpha == 0 ? 0 : -1.0/2.0/beta/beta - 1.0/2.0*d2Alpha*(coshA/sinhA+nBeads/sinhNA)
                + 0.5*dAlpha2*(1/sinhA/sinhA+nBeads*nBeads/sinhNA*coshhNA/sinhNA);

        for (int i=1; i<nBeads; i++) {
            double sinhNmiA = Math.sinh((nBeads-i)*alpha);
            double coshNmiA = Math.cosh((nBeads-i)*alpha);
            double sinhNmip1A = Math.sinh((nBeads-i+1)*alpha);
            double coshNmip1A = Math.cosh((nBeads-i+1)*alpha);

            double sinhRatio = alpha == 0 ? (nBeads-i+1.0)/(nBeads-i) : (sinhNmip1A/sinhNmiA);
            double ki = mass*omegaN*omegaN*sinhRatio;
            chainSigmas[i] = 1.0/Math.sqrt(beta*ki);
            gamma[i] = alpha == 0 ? 1.0/2.0/beta : 1.0/2.0/beta - dAlpha/2.0*(coshNmip1A/sinhNmip1A - (nBeads-i)*sinhA/sinhNmip1A/sinhNmiA);
            dGamma[i] = alpha == 0 ? -1.0/2.0/beta/beta : -1.0/2.0/beta/beta - d2Alpha/2.0*coshNmip1A/sinhNmip1A + (nBeads-i)/2.0*d2Alpha*sinhA/sinhNmip1A/sinhNmiA
                      + dAlpha2/2.0*(nBeads-i+1)/sinhNmip1A/sinhNmip1A
                      + dAlpha2/2.0*(nBeads-i)*coshA/sinhNmip1A/sinhNmiA
                      - dAlpha2/2.0*(nBeads-i)*(nBeads-i+1)*sinhA/sinhNmip1A*coshNmip1A/sinhNmip1A/sinhNmiA
                      - dAlpha2/2.0*(nBeads-i)*(nBeads-i)*sinhA/sinhNmip1A/sinhNmiA*coshNmiA/sinhNmiA;
            f11[i] = alpha == 0 ? (nBeads-i)/(nBeads-i+1.0) : sinhNmiA/sinhNmip1A;
            f1N[i] = alpha == 0 ? 1.0/(nBeads-i+1) : sinhA/sinhNmip1A;

            df11[i] = alpha == 0 ? 0 : dAlpha/sinhNmip1A*((nBeads-i)*coshNmiA-(nBeads-i+1)*sinhNmiA*coshNmip1A/sinhNmip1A);
            df1N[i] = alpha == 0 ? 0 : dAlpha/sinhNmip1A*(coshA - (nBeads-i+1)*sinhA*coshNmip1A/sinhNmip1A);

            d2f11[i] = alpha == 0 ? 0 : d2Alpha/dAlpha*df11[i] + dAlpha2/sinhNmip1A*((nBeads-i)*(nBeads-i)*sinhNmiA
                     - 2*(nBeads-i+1)*(nBeads-i)*coshNmiA*coshNmip1A/sinhNmip1A
                     + 0.5*(nBeads-i+1)*(nBeads-i+1)*sinhNmiA/sinhNmip1A/sinhNmip1A*(3.0+Math.cosh(2*(nBeads-i+1)*alpha)));
            d2f1N[i] = alpha == 0 ? 0 : d2Alpha/dAlpha*df1N[i] + dAlpha2/sinhNmip1A*(sinhA
                     -2*(nBeads-i+1)*coshA*coshNmip1A/sinhNmip1A
                     + 0.5*(nBeads-i+1)*(nBeads-i+1)*sinhA/sinhNmip1A/sinhNmip1A*(3.0 + Math.cosh(2*(nBeads-i+1)*alpha)));
        }
    }

    public double[] getChainSigmas() {
        return chainSigmas;
    }

    public double[] getGamma() {
        return gamma;
    }

    public double[] getdGamma() {
        return dGamma;
    }

    public double[][] getCenterCoefficients() { return new double[][]{f11,f1N}; }

    public double[][] getDCenterCoefficients() {
        return new double[][]{df11,df1N};
    }

    public double[][] getD2CenterCoefficients() {
        return new double[][]{d2f11,d2f1N};
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
        for (int j = 0; j < nBeads; j++) {
            Vector rj = atoms.get(j).getPosition();
            dr.Ev1Mv2(rj, latticePositions[molecule.getIndex()]);
            box.getBoundary().nearestImage(dr);
            uh += 1.0 / nBeads / 2.0 * mass * (omega2 * dr.squared());
            int jj = (j+1)%nBeads;
            Vector rjj = atoms.get(jj).getPosition();
            dr.Ev1Mv2(rjj, rj);
            box.getBoundary().nearestImage(dr);
            uh += 1.0 / 2.0 * mass * (omegaN * omegaN * dr.squared());
        }
        return uh;
    }

    public boolean doTrial() {
        molecule = moleculeSource.getMolecule();

        double uOld = pm.computeOneOldMolecule(molecule);

        IAtomList atoms = molecule.getChildList();
        Vector oldCOM = omega2 == 0 && nGrow == nBeads ? CenterOfMass.position(box, molecule) : null;

        double uhOld = uHarmonic(molecule);
        for (int j = 0; j < nBeads; j++) {
            Vector rj = atoms.get(j).getPosition();
            oldPositions[j].E(rj);
        }

        IAtom atom0 = atoms.get(0);
        Vector prevAtomPosition = atom0.getPosition();
        Vector endAtomPosition = prevAtomPosition;
        int atomStart = 1, numToPlace = nGrow;
        if (nGrow == nBeads) {
            for (int j = 0; j < prevAtomPosition.getD(); j++) {
                prevAtomPosition.setX(j, sigma0 * random.nextGaussian());
            }
            prevAtomPosition.PE(latticePositions[molecule.getIndex()]);
            numToPlace--;
        }
        else {
            int prevAtom = random.nextInt(nBeads);
            atomStart = (prevAtom+1)%nBeads;
            // need to add nBeads because (-1)%nBeads = -1
            prevAtomPosition = atoms.get(prevAtom).getPosition();
            endAtomPosition = atoms.get((atomStart+nGrow)%nBeads).getPosition();
        }
        Vector newCOM = box.getSpace().makeVector();
        for (int k = 0; k < numToPlace; k++) {
            int numPlaced = nBeads - numToPlace + k;
            int kk = (atomStart + k) % nBeads;
            Vector kPosition = atoms.get(kk).getPosition();

            double sigma = chainSigmas[numPlaced];

            for (int j = 0; j < kPosition.getD(); j++) {
                kPosition.setX(j, sigma * random.nextGaussian());
            }

            kPosition.PEa1Tv1( f11[numPlaced], prevAtomPosition);
            kPosition.PEa1Tv1( f1N[numPlaced], endAtomPosition);
            kPosition.PEa1Tv1(1 - f11[numPlaced] - f1N[numPlaced], latticePositions[molecule.getIndex()]);

            prevAtomPosition = kPosition;

            newCOM.PE(kPosition);
        }

        if (nGrow == nBeads && omega2 == 0) {
            newCOM.TE(1.0 / nBeads);
            newCOM.ME(oldCOM);
            for (int k = 0; k < nBeads; k++) {
                atoms.get(k).getPosition().ME(newCOM);
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
