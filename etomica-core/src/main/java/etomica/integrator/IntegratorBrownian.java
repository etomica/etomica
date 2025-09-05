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

import java.io.FileWriter;
import java.io.IOException;

public class IntegratorBrownian extends IntegratorMD {

    protected boolean isLM;
    protected double f2EM;
    protected int numAtoms;
    protected Vector[] rand_0, rand_1;
    protected Vector[] dx;
    protected int n = 0;
    protected FileWriter fw;

    public IntegratorBrownian(PotentialCompute potentialCompute, IRandom random,
                              double h, double temperature, Box box, boolean isLM, double f2EM) {
        super(potentialCompute, random, h, temperature, box);
        this.isLM = isLM;
        this.f2EM = f2EM;
        this.n = 0;
        this.numAtoms = box.getLeafList().size();
        System.out.println(" isLM: " + this.isLM);

        rand_0 = space.makeVectorArray(numAtoms);
        rand_1 = space.makeVectorArray(numAtoms);
        for (int i = 0; i < numAtoms; i++) {
            for (int j = 0; j < space.getD(); j++) {
                rand_0[i].setX(j, random.nextGaussian());
            }
        }
        try {
//            fw = new FileWriter("dx-em");
            fw = new FileWriter("dx-lm");
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    protected void propagatorA(IAtomList atoms, double h) {
        Vector[] forces = potentialCompute.getForces();
        for (int iLeaf = 0; iLeaf < numAtoms; iLeaf++) {
            IAtomKinetic atom = (IAtomKinetic) atoms.get(iLeaf);
            Vector force = forces[iLeaf];
            Vector r = atom.getPosition();

            r.PEa1Tv1(h, force);
//            if (isLM) {
//                r.PEa1Tv1(h, force);
//            } else {
//                r.PEa1Tv1(h/2, force);
//            }
            dx[iLeaf].PEa1Tv1(h , force);
        }
    }

    protected void propagatorO(IAtomList atoms, double h) {
        if (!isothermal) return;
        double sigma = Math.sqrt(2 * temperature * h);
        for (int i = 0; i < numAtoms; i++) {
            for (int j = 0; j < space.getD(); j++) {
                rand_1[i].setX(j, random.nextGaussian());
            }

            IAtomKinetic atom = (IAtomKinetic) atoms.get(i);
            Vector r = atom.getPosition();

            if (isLM) { //LM
                r.PEa1Tv1(sigma / 2, rand_0[i]);
                r.PEa1Tv1(sigma / 2, rand_1[i]);
                dx[i].PEa1Tv1(sigma / 2, rand_0[i]);
                dx[i].PEa1Tv1(sigma / 2, rand_1[i]);
                rand_0[i].E(rand_1[i]);
            } else { //EM
//                sigma = Math.sqrt(2 * temperature * h / 2 * (1 - h * f2EM / 4));
                r.PEa1Tv1(sigma, rand_1[i]);
                dx[i].PEa1Tv1(sigma, rand_1[i]);
            }
        }
    }

    protected void doStepInternal() {
        super.doStepInternal();
        dx = space.makeVectorArray(box.getLeafList().size());
        IAtomList leafList = box.getLeafList();
        propagatorA(leafList, timeStep);
        propagatorO(leafList, timeStep);
        computeForce();

        int n_eq = 100000/5;
//        if (n > n_eq)   System.out.println(n-n_eq-1 + "   " + dx[0].getX(0));
        n++;

//        double sumF2 = 0;
//        Vector[] forces = potentialCompute.getForces();
//        for (int iLeaf = 0; iLeaf < numAtoms; iLeaf++) {
//            sumF2 += forces[iLeaf].squared();
//        }
//        System.out.println(sumF2/(3.0 * numAtoms));

//        double dx2 = 0;
//        for (int i = 0; i < dx.length; i++) {
//            dx2 += dx[i].squared();
//        }
//        dx2 /= (box.getSpace().getD() * box.getLeafList().size());
//        System.out.println(Math.sqrt(dx2));


        int steps = 100000;
        int steps_eq = steps / 5;
        if (n > steps_eq) {
            try {
                fw.write(dx[0].getX(0) + "\n");
//                for (int i=0; i<box.getLeafList().size(); i++) {
//                    fw.write(Math.sqrt(dx[i].squared()) + "\n");
//                }
            }
            catch (IOException ex) {
                throw new RuntimeException(ex);
            }
        }
        if (n == steps_eq + steps) {
            try {
                fw.close();
            }
            catch (IOException ex) {
                throw new RuntimeException(ex);
            }
        }

        dx = space.makeVectorArray(box.getLeafList().size());

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