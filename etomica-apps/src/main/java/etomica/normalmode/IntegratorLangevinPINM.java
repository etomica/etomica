/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorMD;
import etomica.molecule.IMolecule;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.util.Debug;
import etomica.util.random.IRandom;

/**
 * Implements Langevin dynamics using the BAOAB algorithm using path integral normal-mode coordinates.
 * Algorithm originally developed by Leimkuhler and Matthews, Appl. Math. Res. eXpress 2013, 34-56 (2013)
 * and Leimkuhler and Matthews, J. Chem. Phys. 138, 174102 (2013)
 *
 * PI implementation adapted from <a href="https://doi.org/10.1063/1.3489925">Ceriotti et al., J. Chem. Phys. 133,
 * 124104 (2010)</a>, which examined an OBABO integration scheme for PI.
 *
 * Code adapted from <a href="https://github.com/Allen-Tildesley/examples/blob/master/bd_nvt_lj.f90">Allen
 * &amp; Tildesley examples</a>
 *
 * With isothermal = false or gamma = 0, the thermostat will be disabled, resulting in velocity-Verlet
 * integration for NVE.
 *
 * Zeroing net momentum is implemented not by altering the velocities, but by subtracting off the net momentum of mode
 * 0 whenever the transformed coordinates are propagated.  As such, the kinetic energy should be the full 3/2 Nn kT.
 */
public class IntegratorLangevinPINM extends IntegratorMD {

    protected double gamma;
    protected double omega2;
    protected final Vector[] latticePositions;
    protected final double[] mScale;
    protected final MCMoveHO move;

    public IntegratorLangevinPINM(PotentialCompute potentialCompute, IRandom random,
                                  double timeStep, double temperature, Box box, double gamma,
                                  MCMoveHO move, double hbar) {
        super(potentialCompute, random, timeStep, temperature, box);
        setGamma(gamma);
        this.move = move;
        omega2 = this.move.omega2;
        latticePositions = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        for (IMolecule m : box.getMoleculeList()) {
            latticePositions[m.getIndex()].E(m.getChildList().get(0).getPosition());
        }
        int nBeads = box.getMoleculeList().get(0).getChildList().size();
        double mass = box.getLeafList().get(0).getType().getMass();

        mScale = new double[nBeads];

        for (int k = 0; k < nBeads; k++) {
            mScale[k] = this.move.lambda[k]/mass;
        }

        if (omega2 == 0 || nBeads == 1) mScale[nBeads/2] = 1.0;

//        // F = M a;  M = F / a = (sum fi) / (avg ai); fi=1 => sum fi = n
        double[] fq = new double[nBeads];
        for (int k = 0; k < nBeads; k++) {
            for (int i = 0; i < nBeads; i++) {
                fq[k] += this.move.eigenvectors[i][k];
            }
        }

        double[] a = new double[nBeads];
        double mm = box.getLeafList().get(0).getType().getMass();
        double aSum = 0;
        for (int i = 0; i < nBeads; i++) {
            for (int k = 0; k < nBeads; k++) {
                double aq = fq[k] / (mm * mScale[k]);
                a[i] += this.move.eigenvectors[i][k] * aq;
            }
            aSum += a[i];
        }
        // M is the effective mass we have now; we want ring mass = n * atomType mass
        double M = nBeads/(aSum/nBeads);
        double s = (nBeads * box.getLeafList().get(0).getType().getMass()) / M;

        double wn = Math.sqrt(nBeads)*temperature/hbar; // 1/hbar*betan
        System.out.println(" dt-nor-mod = " + Math.sqrt(s));
        System.out.println(" dt-real = " + 1/wn);
        System.out.println();

        for (int i = 0; i < mScale.length; i++) {
            mScale[i] *= s;
        }

        meterKE = new IntegratorPIMD.MeterKineticEnergy(box, mScale);

    }

    public void setGamma(double newGamma) {
        gamma = newGamma;
    }

    protected void propagatorA(double dt) {

        for (IMolecule m : box.getMoleculeList()) {
            IAtomList atoms = m.getChildList();
            Vector[] q = box.getSpace().makeVectorArray(atoms.size());
            int nBeads = atoms.size();
            // first compute collective coordinates
            for (IAtom a : atoms) {
                int i = a.getIndex();
                Vector r = a.getPosition();
                Vector dr = box.getSpace().makeVector();
                dr.Ev1Mv2(r, latticePositions[m.getIndex()]);
                box.getBoundary().nearestImage(dr);
                for (int k = 0; k < nBeads; k++) {
                    q[k].PEa1Tv1(move.eigenvectorsInv[k][i], dr);
                }
            }
            // propagate NM coord

            for (IAtom a : atoms) {
                Vector v = ((IAtomKinetic) a).getVelocity();
                int k = a.getIndex();
                q[k].PEa1Tv1(dt, v);
            }
            // convert to real coord
            for (IAtom a : atoms) {
                int i = a.getIndex();
                Vector r = a.getPosition();
                Vector rOrig = box.getSpace().makeVector();
                rOrig.E(r);
                r.E(latticePositions[m.getIndex()]);
                for (int k = 0; k < nBeads; k++) {
                    r.PEa1Tv1(move.eigenvectors[i][k], q[k]);
                }
                Vector drShift = box.getSpace().makeVector();
                drShift.Ev1Mv2(r, rOrig);
                box.getBoundary().nearestImage(drShift);
                r.Ev1Pv2(rOrig, drShift);
            }
        }
    }

    protected void propagatorB(double dt) {
        Vector[] forces = potentialCompute.getForces();
        for (IMolecule m : box.getMoleculeList()) {
            IAtomList atoms = m.getChildList();
            int nBeads = atoms.size();
            Vector[] fq = box.getSpace().makeVectorArray(nBeads);
            // first compute collective forces
            for (IAtom a : atoms) {
                int i = a.getIndex();
                for (int k = 0; k < nBeads; k++) {
                    fq[k].PEa1Tv1(move.eigenvectors[i][k], forces[a.getLeafIndex()]);
                }
            }

            for (IAtom a : atoms) {
                int i = a.getIndex();
                Vector v = ((IAtomKinetic) a).getVelocity();

                double meff = a.getType().getMass() * mScale[i];
                v.PEa1Tv1(dt / meff, fq[i]);
            }
        }
    }

    protected void propagatorO(double dt) {
        if (!isothermal || gamma == 0) return;
        double c1 = 2.0, c2 = -2.0, c3 = 4.0/3.0, c4 = -2.0/3.0; // Taylor series coefficients

        double x = gamma * dt;
        double c;
        if ( x > 0.0001 ) {
            c = 1-Math.exp(-2.0*x);
        }
        else {
            // Use Taylor expansion for low x
            c = x * (c1 + x * (c2 + x * (c3 + x * c4)));
        }
        c = Math.sqrt ( c );

        double sqrtT = Math.sqrt(temperature);
        Vector rand = space.makeVector();
        IAtomList atoms = box.getLeafList();
        double expX = Math.exp(-x);
        Vector vOld = space.makeVector();
        for (IAtom atom : atoms) {
            IAtomKinetic a = (IAtomKinetic) atom;
            double m = a.getType().getMass() * mScale[a.getIndex()];
            double sqrtM = Math.sqrt(m);
            for (int i = 0; i < rand.getD(); i++) {
                rand.setX(i, c * sqrtT / sqrtM * random.nextGaussian());
            }
            Vector v = a.getVelocity();
            vOld.E(v);
            v.TE(expX);
            v.PE(rand);
        }
    }

    protected void doStepInternal() {
        super.doStepInternal();

        propagatorB(timeStep/2);
        propagatorA(timeStep/2);
        propagatorO(timeStep);
        propagatorA(timeStep/2);
        computeForce();
        propagatorB(timeStep/2);

        computeKE();
    }

    protected void computeKE() {
        IAtomList atoms = box.getLeafList();
        currentKineticEnergy = 0;

        for (IAtom atom : atoms) {
            IAtomKinetic a = (IAtomKinetic) atom;
            Vector velocity = a.getVelocity();
            currentKineticEnergy += 0.5 * a.getType().getMass() * mScale[a.getIndex()] * velocity.squared();
        }
    }

    public void reset() {
        super.reset();

        if (Debug.ON && Debug.DEBUG_NOW) {
            IAtomList pair = Debug.getAtoms(box);
            if (pair != null) {
                Vector dr = space.makeVector();
                dr.Ev1Mv2(pair.get(1).getPosition(), pair.get(0).getPosition());
                System.out.println(pair + " dr " + dr);
            }
        }

        computeForce();
    }

    public void computeForce() {
        eventManager.forcePrecomputed();

        currentPotentialEnergy = potentialCompute.computeAll(true);
        eventManager.forceComputed();
    }

    public void randomizeMomenta() {
        if (thermostatNoDrift) {
            throw new RuntimeException("Langevin for PI is currently unable to prevent drift");
        }
        atomActionRandomizeVelocity.setTemperature(temperature);
        IAtomList leafList = box.getLeafList();
        for (IAtom a : leafList) {
            atomActionRandomizeVelocity.actionPerformed(a);
            ((IAtomKinetic) a).getVelocity().TE(1 / Math.sqrt(mScale[a.getIndex()]));
        }
        if (alwaysScaleMomenta) {
            scaleMomenta();
        }
    }

    public void shiftMomenta() {
        // do nothing
    }

    protected void scaleMomenta(Vector t) {
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.size();
        currentKineticEnergy = 0;
        if (nLeaf == 0) return;
        // calculate current kinetic temperature.
        for (int i = 0; i < space.D(); i++) {
            // scale independently in each dimension
            double sum = 0.0;
            int nLeafNotFixed = 0;
            for (IAtom value : leafList) {
                IAtomKinetic atom = (IAtomKinetic) value;
                double mass = atom.getType().getMass();
                if (mass == Double.POSITIVE_INFINITY) continue;
                double v = atom.getVelocity().getX(i);
                sum += mass * mScale[atom.getIndex()] * v * v;
                nLeafNotFixed++;
            }
            if (sum == 0) {
                if (t.getX(i) == 0) {
                    continue;
                }
                if (i > 0) {
                    // wonky.  possible in theory.  but then, you called
                    // scaleMomenta, so you're probably a bad person and
                    // deserve this.
                    throw new RuntimeException("atoms have no velocity component in " + i + " dimension");
                }
                randomizeMomenta();
                i--;
                // try again, we could infinite loop in theory
                continue;
            }
            double s = Math.sqrt(t.getX(i) / (sum / nLeafNotFixed));
            currentKineticEnergy += 0.5 * sum * s * s;
            if (s == 1) continue;
            for (IAtom value : leafList) {
                IAtomKinetic atom = (IAtomKinetic) value;
                Vector vel = atom.getVelocity();
                vel.setX(i, vel.getX(i) * s); //scale momentum
            }
        }
    }


    public void postRestore() {
        super.postRestore();
        computeForce();
    }
}