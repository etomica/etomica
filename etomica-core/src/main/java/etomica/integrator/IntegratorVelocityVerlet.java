/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.atom.AtomSetSinglet;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.util.Debug;
import etomica.util.random.IRandom;

public class IntegratorVelocityVerlet extends IntegratorMD {

    public IntegratorVelocityVerlet(PotentialCompute potentialCompute, IRandom random,
                                    double timeStep, double temperature, Box box) {
        super(potentialCompute, random, timeStep, temperature, box);
    }

//--------------------------------------------------------------
// steps all particles across time interval tStep

    // assumes one box
    protected void doStepInternal() {
        super.doStepInternal();
        if (Debug.ON && Debug.DEBUG_NOW) {
            IAtomList pair = Debug.getAtoms(box);
            if (pair != null) {
                Vector dr = space.makeVector();
                dr.Ev1Mv2(pair.get(1).getPosition(), pair.get(0).getPosition());
                System.out.println(pair + " dr " + dr);
            }
        }
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.size();
        Vector[] forces = potentialCompute.getForces();
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic) leafList.get(iLeaf);
            Vector force = forces[iLeaf];
            Vector r = a.getPosition();
            Vector v = a.getVelocity();
            if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                System.out.println("first " + a + " r=" + r + ", v=" + v + ", f=" + force);
            }
            v.PEa1Tv1(0.5 * timeStep * a.getType().rm(), force);  // p += f(old)*dt/2
            r.PEa1Tv1(0.5 * timeStep, v);         // r += p*dt/m/2
        }

        if (isothermal) {
            doThermostatInternal();
        }

        //Finish integration step
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic) leafList.get(iLeaf);
            Vector r = a.getPosition();
            Vector v = a.getVelocity();
            if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                System.out.println("second " + a + " v=" + v + ", f=" + forces[iLeaf]);
            }
            r.PEa1Tv1(0.5 * timeStep, v);         // r += p*dt/m/2
        }

        eventManager.forcePrecomputed();

        currentPotentialEnergy = potentialCompute.computeAll(true);

        eventManager.forceComputed();

        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic) leafList.get(iLeaf);
            Vector v = a.getVelocity();
            if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                System.out.println("third " + a + " v=" + v + ", f=" + forces[iLeaf]);
            }
            v.PEa1Tv1(0.5 * timeStep * a.getType().rm(), forces[iLeaf]);  //p += f(new)*dt/2
        }

        eventManager.finalizeStep();

        currentKineticEnergy = 0;
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic) leafList.get(iLeaf);
            currentKineticEnergy += 0.5 * a.getType().getMass() * a.getVelocity().squared();
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

    public void postRestore() {
        super.postRestore();
        precomputeForce();
    }
}
