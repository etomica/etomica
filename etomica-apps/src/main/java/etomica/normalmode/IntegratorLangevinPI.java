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
 * Implements Langevin dynamics using the BAOAB algorithm using path integral staging coordinates.
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
public class IntegratorLangevinPI extends IntegratorMD {

    protected double gamma;
    protected final double[] sigma, f11, f1N;
    protected double omega2HO;
    protected final Vector[] latticePositions;
    protected final double[] mScale, fScale, fScale0;

    public IntegratorLangevinPI(PotentialCompute potentialCompute, IRandom random,
                                double timeStep, double temperature, Box box, double gamma,
                                MCMoveHOReal2 move, double hbar, double omega2HO) {
        super(potentialCompute, random, timeStep, temperature, box);
        setGamma(gamma);

        sigma = move.chainSigmas;
        f11 = move.f11;
        f1N = move.f1N;
        this.omega2HO = omega2HO;
        latticePositions = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        for (IMolecule m : box.getMoleculeList()) {
            latticePositions[m.getIndex()].E(m.getChildList().get(0).getPosition());
        }
        int nBeads = box.getMoleculeList().get(0).getChildList().size();
        double omegaN = Math.sqrt(nBeads)/(hbar* move.beta);
        double omegaN2 = omegaN*omegaN;
        double D = 2 + move.omega2 / (nBeads*omegaN2);
        double alpha = Math.log(D/2 + Math.sqrt(D*D/4 - 1));
        mScale = new double[nBeads];
        fScale = new double[nBeads];
        fScale0 = new double[nBeads];

        mScale[0] = (alpha == 0 || nBeads == 1) ? nBeads*omegaN2/omega2HO : 2.0*nBeads*omegaN2/omega2HO*Math.sinh(alpha)*Math.tanh(nBeads*alpha/2.0);

// Check how m0 varies with x for the EC-stage
//        System.out.println(mScale[0]/nBeads);
//        if (alpha != 0) {
//            for (int i=1;i<=100;i++) {
//                double x = 1.0*i;
//                double n = 20 * x;
//                double d = 1 + 0.5 * (x / n) * (x / n);
//                double a = Math.log(d + Math.sqrt(d * d - 1));
//                double m0 = 1 * (2 / n) * (n / x) * (n / x) * Math.sinh(a) * Math.tanh(n * a / 2);
//                System.out.println(x + "  " + m0);
//            }
//            System.exit(0);
//        }

        fScale0[0] = 1;
        for (int i=1; i<nBeads; i++) {
            fScale0[i] = alpha == 0 ? 1.0 : Math.cosh((nBeads / 2.0 - i)*alpha) / Math.cosh(nBeads/2.0*alpha);
            fScale[i]  = alpha == 0 ? (nBeads - i - 1.0)/(nBeads - i) : Math.sinh((nBeads - i - 1) * alpha)/Math.sinh((nBeads - i)*alpha);
            mScale[i]  = alpha == 0 ? nBeads*omegaN2/omega2HO*(nBeads - i + 1.0)/(nBeads - i) : nBeads*omegaN2/omega2HO*Math.sinh((nBeads - i + 1) * alpha)/Math.sinh((nBeads - i)*alpha);
        }

      // The "s" rescaling is no longer needed as s=1 always!
        // F = M a;  M = F / a = (sum fi) / (avg ai); fi=1 => sum fi = nBeads
//        double[] fu = new double[nBeads];
//        for (int i=nBeads-1; i>=0; i--) {
//            fu[0] += fScale0[i];
//            if (i>0) {
//                fu[i] = 1;
//                if (i < nBeads - 1) {
//                    fu[i] += fScale[i] * fu[i + 1];
//                }
//            }
//        }
//
//        double[] a = new double[nBeads];
//        double mm = box.getLeafList().get(0).getType().getMass();
//        double aSum = 0;
//        for (int i=0; i<nBeads; i++) {
//            a[i] = fu[i] / (mm*mScale[i]);
//            if (i>0) {
//                a[i] += f11[i] * a[i-1];
//                a[i] += f1N[i] * a[0];
//            }
//            aSum += a[i];
//        }
//        // M is the effective mass we have now; we want ring mass = nBeads * atomType mass
//        double M = nBeads/(aSum/nBeads);
//        double s = (nBeads * box.getLeafList().get(0).getType().getMass()) / M;


        // mi need to be larger to match the real COM oscillations
        if (alpha == 0) {
            double s = omega2HO/omegaN/omegaN*(1 + 1.0/12.0*(nBeads*nBeads-1)/nBeads);
            for (int i = 0; i < mScale.length; i++) {
                mScale[i] *= s;
            }
        }

        meterKE = new IntegratorPIMD.MeterKineticEnergy(box, mScale);
    }

    public void setGamma(double newGamma) {
        gamma = newGamma;
    }

    protected void propagatorA(double dt) {
        Vector drPrev0 = box.getSpace().makeVector();
        Vector drPrev = box.getSpace().makeVector();

        for (IMolecule m : box.getMoleculeList()) {
            IAtomList atoms = m.getChildList();
            Vector dr0 = box.getSpace().makeVector();
            dr0.Ev1Mv2(atoms.get(0).getPosition(), latticePositions[m.getIndex()]);
            box.getBoundary().nearestImage(dr0);
            Vector[] u = box.getSpace().makeVectorArray(atoms.size());
            // first compute collective coordinates
            for (IAtom a : atoms) {
                int i = a.getIndex();
                Vector r = a.getPosition();
                u[i].Ev1Mv2(r, latticePositions[m.getIndex()]);
                box.getBoundary().nearestImage(u[i]);
                Vector usave = box.getSpace().makeVector();
                usave.E(u[i]);
                if (i > 0) {
                    u[i].PEa1Tv1(-f11[i], drPrev0);
                    u[i].PEa1Tv1(-f1N[i], dr0);
                }
                drPrev0.E(usave);

                Vector v = ((IAtomKinetic)a).getVelocity();
                u[i].PEa1Tv1(dt, v);

                Vector rOrig = box.getSpace().makeVector();
                rOrig.E(r);
                r.E(latticePositions[m.getIndex()]);
                if (i>0) {
                    r.PEa1Tv1(f11[i], drPrev);
                    r.PEa1Tv1(f1N[i], u[0]);
                }
                r.PE(u[i]);
                drPrev.Ev1Mv2(r, latticePositions[m.getIndex()]);
                // shift atom back to side of box where it started, even if outside
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
            int n = atoms.size();
            Vector[] fu = box.getSpace().makeVectorArray(n);

            // first compute collective forces
            for (int i = n - 1; i >= 0; i--) {
                IAtom a = atoms.get(i);
                Vector foo = box.getSpace().makeVector();
                foo.E(forces[a.getLeafIndex()]);
                fu[0].PEa1Tv1(fScale0[i], foo);
                if (i > 0) {
                    fu[i].E(foo);
                    if (i < n - 1) {
                        fu[i].PEa1Tv1(fScale[i], fu[i + 1]);
                    }
                }
            }
            for (IAtom a : atoms) {
                int i = a.getIndex();
                Vector v = ((IAtomKinetic) a).getVelocity();

                double meff = a.getType().getMass() * mScale[i];
                v.PEa1Tv1(dt / meff, fu[i]);
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
