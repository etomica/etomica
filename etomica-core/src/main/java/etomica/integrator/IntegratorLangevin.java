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

    public IntegratorLangevin(PotentialCompute potentialCompute, IRandom random,
                              double timeStep, double temperature, Box box, double gamma) {
        super(potentialCompute, random, timeStep, temperature, box);
        setGamma(gamma);
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
        Vector netMomentum = space.makeVector();
        double totalMass = 0;
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
                totalMass += m;
                netMomentum.PEa1Tv1(m, rand);
            }
        }
        if (thermostatNoDrift) {
            netMomentum.TE(1.0/totalMass);
            for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
                IAtomKinetic a = (IAtomKinetic) atoms.get(iLeaf);
                a.getVelocity().ME(netMomentum);
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
        atomActionRandomizeVelocity.setTemperature(temperature);
        IAtomList leafList = box.getLeafList();
        for (IAtom a : leafList) {
            atomActionRandomizeVelocity.actionPerformed(a);
            ((IAtomKinetic) a).getVelocity().TE(1 / Math.sqrt(1.0));
        }
        if (alwaysScaleMomenta) {
            scaleMomenta();
        }
    }

    public void postRestore() {
        super.postRestore();
        computeForce();
    }
}