/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.util.Debug;
import etomica.util.random.IRandom;

/**
 * Implements Langevin dynamics using the BAOAB algorithm.
 * originally developed by Leimkuhler and Matthews, Appl. Math. Res. eXpress 2013, 34–56 (2013)
 * and Leimkuhler and Matthews, J. Chem. Phys. 138, 174102 (2013)
 *
 * code adapted from <a href="https://github.com/Allen-Tildesley/examples/blob/master/bd_nvt_lj.f90">...</a>
 *
 * With isothermal = false or gamma = 0, the integration thermostat will be disabled, resulting
 * in velocity-Verlet integration for NVE.
 */
public class IntegratorLangevin extends IntegratorMD {

    protected double gamma;
    protected final Vector netMomentum;
    protected double totalMass0;


    public IntegratorLangevin(PotentialCompute potentialCompute, IRandom random,
                              double timeStep, double temperature, Box box, double gamma) {
        super(potentialCompute, random, timeStep, temperature, box);
        setGamma(gamma);
        netMomentum = space.makeVector();
    }

    public void setGamma(double newGamma) {
        gamma = newGamma;
    }

    protected void propagatorA(IAtomList atoms, double dt) {
        int nLeaf = atoms.size();
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic) atoms.get(iLeaf);
            Vector r = a.getPosition();
            Vector v = a.getVelocity();
            r.PEa1Tv1(dt, v);

            if (thermostatNoDrift) v.PEa1Tv1(-dt/totalMass0, netMomentum);

        }
    }

    protected void propagatorB(IAtomList atoms, double dt) {
        int nLeaf = atoms.size();
        Vector[] forces = potentialCompute.getForces();
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic) atoms.get(iLeaf);
            Vector force = forces[iLeaf];
            Vector v = a.getVelocity();
            v.PEa1Tv1(dt * a.getType().rm(), force);  // p += f(old)*dt/2
            if (thermostatNoDrift) {
                netMomentum.PEa1Tv1(dt, force);
            }

        }
    }

    protected void propagatorO(IAtomList atoms, double dt) {
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
        int nLeaf = atoms.size();
        double expX = Math.exp(-x);
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic) atoms.get(iLeaf);
            double m = a.getType().getMass();
            double sqrtM = Math.sqrt(m);
            for (int i=0; i<rand.getD(); i++) {
                rand.setX(i, c * sqrtT / sqrtM * random.nextGaussian());
            }
            Vector v = a.getVelocity();
            v.TE(expX);
            v.PE(rand);



            if (thermostatNoDrift) {
                netMomentum.PEa1Tv1(-m, v);
            }

        }
    }



    protected void doStepInternal() {
        super.doStepInternal();

        IAtomList leafList = box.getLeafList();
        propagatorB(leafList, timeStep/2);
        propagatorA(leafList, timeStep/2);
        propagatorO(leafList, timeStep);
        propagatorA(leafList, timeStep/2);
        computeForce();
        propagatorB(leafList, timeStep/2);

        computeKE(leafList);
    }

    protected void computeKE(IAtomList atoms) {
        currentKineticEnergy = 0;

        int nLeaf = atoms.size();
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic) atoms.get(iLeaf);
            Vector velocity = a.getVelocity();
            currentKineticEnergy += 0.5 * a.getType().getMass() * velocity.squared();
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
//        if (thermostatNoDrift) {
//            throw new RuntimeException("Langevin for PI is currently unable to prevent drift");
//        }
        atomActionRandomizeVelocity.setTemperature(temperature);
        IAtomList leafList = box.getLeafList();
        for (IAtom a : leafList) {
            atomActionRandomizeVelocity.actionPerformed(a);
//            ((IAtomKinetic) a).getVelocity().TE(1 / Math.sqrt(mScale[a.getIndex()]));
            ((IAtomKinetic) a).getVelocity().TE(1 / Math.sqrt(1.0));
        }
        if (alwaysScaleMomenta) {
            scaleMomenta();
        }
    }

    public void shiftMomenta() {
        netMomentum.E(0);
        totalMass0 = 0;
        for (IAtom a : box.getLeafList()) {
            int i = a.getIndex();
            if (i != 0) continue;
            double m = a.getType().getMass();
            netMomentum.PEa1Tv1(m, ((IAtomKinetic) a).getVelocity());
            totalMass0 += m;
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
            for (IAtom value : leafList) {
                IAtomKinetic atom = (IAtomKinetic) value;
                double mass = atom.getType().getMass();
                if (mass == Double.POSITIVE_INFINITY) continue;
                double v = atom.getVelocity().getX(i);
                sum += mass * v * v;
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