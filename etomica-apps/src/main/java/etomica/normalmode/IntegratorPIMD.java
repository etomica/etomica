/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorMD;
import etomica.molecule.IMolecule;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.units.dimensions.Energy;
import etomica.util.Debug;
import etomica.util.random.IRandom;

public class IntegratorPIMD extends IntegratorMD {

    protected final double[] sigma, f11, f1N;
    protected double omega2;
    protected final Vector[] latticePositions;
    protected final double[] mScale, fScale, fScale0;

    /**
     * Constructs integrator with a default for non-isothermal sampling.
     *
     * @param potentialCompute PotentialMaster instance used to compute the energy and forces
     * @param random           random number generator used for initial velocities and some thermostats
     * @param timeStep         time step for integration
     * @param temperature      used by thermostat and/or to initialize velocities
     * @param box
     */
    public IntegratorPIMD(PotentialCompute potentialCompute, IRandom random, double timeStep, double temperature, Box box, MCMoveHOReal2 move, double hbar) {
        super(potentialCompute, random, timeStep, temperature, box);
        sigma = move.chainSigmas;
        f11 = move.f11;
        f1N = move.f1N;
        omega2 = move.omega2;
        latticePositions = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        for (IMolecule m : box.getMoleculeList()) {
            latticePositions[m.getIndex()].E(m.getChildList().get(0).getPosition());
        }
        int n = box.getMoleculeList().get(0).getChildList().size();
        double omegaN = temperature*n/hbar;
        double D = 2 + omega2 / (omegaN*omegaN);
        double alpha = Math.log(D/2 + Math.sqrt(D*D/4 - 1));
        mScale = new double[n];
        fScale = new double[n];
        fScale0 = new double[n];
        mScale[0] = 2.0*Math.sinh(alpha) * Math.tanh(n*alpha/2.0);
//        mScale[0] = 2.0*Math.sinh(alpha) * Math.tanh(n*alpha/2.0)*omegaN*omegaN/omega2;
        if (alpha == 0 || n == 1) mScale[0] = 1.0;
        fScale0[0] = 1;
        for (int i=1; i<n; i++) {
            fScale0[i] = alpha == 0 ? 1.0 : Math.cosh((n / 2.0 - i)*alpha) / Math.cosh(n/2.0*alpha);
            fScale[i]  = alpha == 0 ? (n - i - 1.0)/(n - i) : (Math.sinh((n - i - 1) * alpha) / Math.sinh((n - i)*alpha));
            mScale[i]  = alpha == 0 ? (n - i + 1.0)/(n - i) : (Math.sinh((n - i + 1) * alpha) / Math.sinh((n - i)*alpha));
//            mScale[i]  = alpha == 0 ? (n - i + 1.0)/(n - i) : (Math.sinh((n - i + 1) * alpha) / Math.sinh((n - i)*alpha)*omegaN*omegaN/omega2);
        }

        // F = M a;  M = F / a = (sum fi) / (sum ai); fi=1 => sum fi = n
        double[] fu = new double[n];
        for (int i=n-1; i>=0; i--) {
            fu[0] += fScale0[i];
            if (i>0) {
                fu[i] = 1;
                if (i < n - 1) {
                    fu[i] += fScale[i] * fu[i + 1];
                }
            }
        }

        double[] ai = new double[n];
        double mm = box.getLeafList().get(0).getType().getMass();
        double aSum = 0;
        for (int i=0; i<n; i++) {
            ai[i] = fu[i] / (mm*mScale[i]);
            if (i>0) {
                ai[i] += f11[i] * ai[i-1];
                ai[i] += f1N[i] * ai[0];
            }
            aSum += ai[i];
        }
        // M is the effective mass we have now; we want ring mass = n * atomType mass
        double M = n / (aSum/n);
        double s = (n * box.getLeafList().get(0).getType().getMass()) / M;
        for (int i=0; i<mScale.length; i++) {
            mScale[i] *= s;
        }

        meterKE = new MeterKineticEnergy(box, mScale);
    }

    public void doStepInternal() {
        super.doStepInternal();
        Vector[] forces = potentialCompute.getForces();
        int n = box.getMoleculeList().get(0).getChildList().size();
        Vector drPrev0 = box.getSpace().makeVector();
        Vector drPrev = box.getSpace().makeVector();

        for (IMolecule m : box.getMoleculeList()) {
            IAtomList atoms = m.getChildList();
            Vector dr0 = box.getSpace().makeVector();
            dr0.Ev1Mv2(atoms.get(0).getPosition(), latticePositions[m.getIndex()]);
            box.getBoundary().nearestImage(dr0);
            Vector[] u = box.getSpace().makeVectorArray(atoms.size());
            Vector[] fu = box.getSpace().makeVectorArray(atoms.size());
            // first compute collective coordinates and forces
            for (int i=atoms.size()-1; i>=0; i--) {
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
            // now propogate coordinates and velocities.  collective velocities are stored
            // as atom's velocity.
            for (IAtom a : atoms) {
                // velocity-Verlet: propagate velocity by half, u by full timestep
//                v.PEa1Tv1(0.5 * timeStep * a.getType().rm(), force);
//                r.PEa1Tv1(timeStep, v);
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
                double meff = a.getType().getMass() * mScale[i];
                v.PEa1Tv1(0.5 * timeStep / meff, fu[i]);
                u[i].PEa1Tv1(timeStep, v);

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

        eventManager.forcePrecomputed();

        currentPotentialEnergy = potentialCompute.computeAll(true);

        eventManager.forceComputed();

        //Finish integration step
        for (IMolecule m : box.getMoleculeList()) {
            IAtomList atoms = m.getChildList();
            double mass = atoms.get(0).getType().getMass();
            // first compute collective forces
            Vector[] fu = box.getSpace().makeVectorArray(atoms.size());
            for (int i=atoms.size()-1; i>=0; i--) {
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
                // velocity-Verlet: propagate velocity by half, u by full timestep
//                velocity.PEa1Tv1(0.5 * timeStep * a.getType().rm(), forces[iLeaf]);
                int i = a.getIndex();
                Vector v = ((IAtomKinetic) a).getVelocity();
                double meff = mass * mScale[i];
                v.PEa1Tv1(0.5 * timeStep / meff, fu[i]);
            }
        }

        eventManager.finalizeStep();

        currentKineticEnergy = 0;
        for (IMolecule m : box.getMoleculeList()) {
            IAtomList atoms = m.getChildList();
            double mass = atoms.get(0).getType().getMass();
            for (int i = 0; i < atoms.size(); i++) {
                double meff = mass * mScale[i];
                IAtom a = atoms.get(i);
                Vector velocity = ((IAtomKinetic) a).getVelocity();
                currentKineticEnergy += 0.5 * meff * velocity.squared();
            }
        }

        if (isothermal) {
            doThermostatInternal();
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

        precomputeForce();
    }

    public void precomputeForce() {
        eventManager.forcePrecomputed();

        currentPotentialEnergy = potentialCompute.computeAll(true);
        eventManager.forceComputed();
    }

    public double[] getMassScale() {
        return mScale;
    }

    public void randomizeMomenta() {
        atomActionRandomizeVelocity.setTemperature(temperature);
        IAtomList leafList = box.getLeafList();
        for (IAtom a : leafList) {
            atomActionRandomizeVelocity.actionPerformed(a);
            ((IAtomKinetic)a).getVelocity().TE(1/Math.sqrt(mScale[a.getIndex()]));
        }
        if (alwaysScaleMomenta) {
            if (thermostatNoDrift) {
                shiftMomenta();
            }
            scaleMomenta();
        }
    }

    public void shiftMomenta() {
        Vector[] p = box.getSpace().makeVectorArray(box.getMoleculeList().get(0).getChildList().size());
        double[] totalMass = new double[p.length];
        for (IAtom a : box.getLeafList()) {
            int i = a.getIndex();
            p[i].PEa1Tv1(a.getType().getMass(), ((IAtomKinetic)a).getVelocity());
            totalMass[i] += a.getType().getMass();
        }
        for (int i=0; i<p.length; i++) {
            p[i].TE(1/totalMass[i]);
        }
        for (IAtom a : box.getLeafList()) {
            ((IAtomKinetic)a).getVelocity().ME(p[a.getIndex()]);
        }
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
            for (int iAtom = 0; iAtom < nLeaf; iAtom++) {
                IAtomKinetic atom = (IAtomKinetic) leafList.get(iAtom);
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
            for (int iAtom = 0; iAtom < nLeaf; iAtom++) {
                IAtomKinetic atom = (IAtomKinetic) leafList.get(iAtom);
                Vector vel = atom.getVelocity();
                vel.setX(i, vel.getX(i) * s); //scale momentum
            }
        }
    }


    public static class MeterKineticEnergy extends DataSourceScalar {

        protected final Box box;
        protected final double[] mScale;

        public MeterKineticEnergy(Box box, double[] mScale) {
            super("KE", Energy.DIMENSION);
            this.box = box;
            this.mScale = mScale;
        }

        public double getDataAsScalar() {
            double ke = 0.0;
            for (IAtom atom : box.getLeafList()) {
                double mass = atom.getType().getMass() * mScale[atom.getIndex()];
                if (mass == Double.POSITIVE_INFINITY) continue;
                ke += 0.5 * mass * (((IAtomKinetic)atom).getVelocity().squared());
            }
            return ke;
        }
    }
}
