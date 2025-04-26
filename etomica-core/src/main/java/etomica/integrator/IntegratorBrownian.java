/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.util.Debug;
import etomica.util.random.IRandom;

public class IntegratorBrownian extends IntegratorMD {

    protected boolean isLM;
    public IntegratorBrownian(PotentialCompute potentialCompute, IRandom random,
                              double h, double temperature, Box box) {
        super(potentialCompute, random, h, temperature, box);
        this.isLM = true;
        System.out.println(" isLM: " + isLM);
    }


    protected void propagatorA(IAtomList atoms, double h) {
        int nLeaf = atoms.size();
        Vector[] forces = potentialCompute.getForces();
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic) atoms.get(iLeaf);
            Vector force = forces[iLeaf];
            Vector r = a.getPosition();
            r.PEa1Tv1(h/a.getType().getMass(), force);
        }
    }

    protected void propagatorO(IAtomList atoms, double h) {
        if (!isothermal) return;
        Vector rand = space.makeVector();
        int nLeaf = atoms.size();
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic) atoms.get(iLeaf);
            double mass = a.getType().getMass();
            double sqrtT = Math.sqrt(2*temperature*h/mass);
            for (int i=0; i<rand.getD(); i++) {
                rand.setX(i, random.nextGaussian());
            }
            Vector r = a.getPosition();
            if (isLM) { //LM
                r.PEa1Tv1(sqrtT/2, rand);
                r.PEa1Tv1(sqrtT/2*Math.sqrt(mass/temperature), a.getVelocity());
                a.getVelocity().Ea1Tv1(Math.sqrt(temperature/mass), rand);
            } else { //EM
                r.PEa1Tv1(sqrtT, rand);
            }
        }
    }

    protected void doStepInternal() {
        super.doStepInternal();

        IAtomList leafList = box.getLeafList();
        propagatorA(leafList, timeStep);
        propagatorO(leafList, timeStep);
        computeForce();

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

    public void postRestore() {
        super.postRestore();
        computeForce();
    }
}