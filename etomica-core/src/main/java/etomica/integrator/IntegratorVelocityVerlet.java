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
import etomica.space3d.Vector3D;
import etomica.util.Debug;
import etomica.util.random.IRandom;

import java.util.ArrayList;

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
        double dotValF =0;
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
        ArrayList<double[]> rFinal = new ArrayList<>();
        ArrayList<double[]> rInitial = new ArrayList<>();
        ArrayList<double[]> fFinal = new ArrayList<>();
        ArrayList<double[]> rDiff = new ArrayList<>();
        ArrayList<Double> rdotF = new ArrayList<>();
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic) leafList.get(iLeaf);
            Vector force = forces[iLeaf];
            Vector r = a.getPosition();
            Vector v = a.getVelocity();
            if (Debug.ON && Debug.DEBUG_NOW ) {
                rInitial.add(r.toArray());
                System.out.println("first " + a + " r=" + r + ", v=" + v + ", f=" + force);
            }
            Vector rOld = a.getPosition();
            v.PEa1Tv1(0.5 * timeStep * a.getType().rm(), force);  // p += f(old)*dt/2
            r.PEa1Tv1(timeStep, v);         // r += p*dt/m
         /*  System.out.println("rdiff  " + rdiff);
           rdiff.dot(force);
            System.out.println("prod " +rdiff );*/
        }

        eventManager.forcePrecomputed();

        currentPotentialEnergy = potentialCompute.computeAll(true);



        eventManager.forceComputed();

        //Finish integration step
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic) leafList.get(iLeaf);
            Vector velocity = a.getVelocity();
            Vector r = a.getPosition();
            if (Debug.ON  && Debug.DEBUG_NOW ) {
                rFinal.add(r.toArray());
                fFinal.add(forces[iLeaf].toArray());
                System.out.println("second " + a +" r=" + r + ", v=" + velocity + ", f=" + forces[iLeaf]);
            }
            velocity.PEa1Tv1(0.5 * timeStep * a.getType().rm(), forces[iLeaf]);  //p += f(new)*dt/2
        }
        eventManager.preThermostat();
        currentKineticEnergy = 0;
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic) leafList.get(iLeaf);
            Vector velocity = a.getVelocity();
            currentKineticEnergy += 0.5 * a.getType().getMass() * velocity.squared();
        }

        if (isothermal) {
            doThermostatInternal();
        }
      /*  IAtomKinetic a = (IAtomKinetic) leafList.get(0);
        for (int i=0; i<fFinal.size(); i++){
            Vector v1 = Vector.of(rInitial.get(i));
            Vector v2 = Vector.of(rFinal.get(i));
            v1.ME(v2);
            rDiff.add(v1.toArray());
            Vector f1 = Vector.of(fFinal.get(i));
            double dotVal = v1.dot(f1);
            dotValF += dotVal;
            rdotF.add(dotVal);
        }*/
      // System.out.println("Step: " + stepCount+" PE: "+ currentPotentialEnergy + " KE: " + currentKineticEnergy + " TE: "+ (currentKineticEnergy+currentPotentialEnergy)  +" "+ dotValF+ "\n");
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
