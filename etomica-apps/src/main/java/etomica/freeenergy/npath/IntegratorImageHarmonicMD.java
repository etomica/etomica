/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.freeenergy.npath;

import etomica.atom.AtomSetSinglet;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.chem.elements.Iron;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.PotentialCalculationForcePressureSum;
import etomica.potential.PotentialMaster;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.Debug;
import etomica.util.random.IRandom;

/**
 * Created by andrew on 4/30/17.
 */
public class IntegratorImageHarmonicMD extends IntegratorVelocityVerlet {

    protected final Vector dr, drTmp, dv, dvTmp, df;
    protected final Vector fTot, vTot;
    protected P1ImageHarmonic p1;
    protected Boundary boundary;

    public IntegratorImageHarmonicMD(PotentialMaster potentialMaster, IRandom random, double timeStep, double temperature, Space space) {
        super(potentialMaster, random, timeStep, temperature, space);
        dr = space.makeVector();
        drTmp = space.makeVector();
        dv = space.makeVector();
        dvTmp = space.makeVector();
        df = space.makeVector();
        fTot = space.makeVector();
        vTot = space.makeVector();
    }

    public void setP1Harmonic(P1ImageHarmonic p1) {
        this.p1 = p1;
    }

    public void setBox(Box box) {
        super.setBox(box);
        boundary = box.getBoundary();
    }

    public double getRandomizeProbability() {
        double w = p1.getW();
        double rm2 = 2 / Iron.INSTANCE.getMass();
        double o = Math.sqrt(2 * w * rm2);
        double x = o * timeStep / Math.PI;
        return 10 * x * x;
    }

    protected void doStepInternal() {
        // from IntegratorMD
        currentTime += timeStep;
        // IntegratorVelocityVerlet code
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            MyAgent agent = agentManager.getAgent(a);
            Vector r = a.getPosition();
            Vector v = a.getVelocity();
            if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                System.out.println("first "+a+" r="+r+", v="+v+", f="+agent.force);
            }
            v.PEa1Tv1(0.5 * timeStep * a.getType().rm(), agent.force);  // p += f(old)*dt/2
        }

        double omega = 0;
        double sinomegat = 0;
        double cosomegat = 0;
        Vector offset = p1.getOffset();
        double pRandomize = getRandomizeProbability();
        // if w is very strong, our springs will oscillate multiple times
        // during our timestep, which is mostly useless.  If the frequency
        // happens to match our timestep, then the dimer internal coordinates
        // will never be sampled.  Break this by randomizing
        boolean doRandomize = pRandomize > 1 || pRandomize > random.nextDouble();
        int passes = doRandomize ? 2 : 1;
        double tStep0 = doRandomize ? timeStep * 0.5 : timeStep;
        double w = p1.getW();
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            int iLeaf1 = p1.getPartner(iLeaf);
            if (iLeaf1 < iLeaf) {
                continue;
            }
            IAtomKinetic atom0 = (IAtomKinetic) leafList.getAtom(iLeaf);
            IAtomKinetic atom1 = (IAtomKinetic) leafList.getAtom(iLeaf1);
            Vector r0 = atom0.getPosition();
            Vector r1 = atom1.getPosition();
            Vector v0 = atom0.getVelocity();
            Vector v1 = atom1.getVelocity();

            vTot.Ev1Pv2(v1, v0);
            vTot.TE(0.5);
            r0.PEa1Tv1(timeStep, vTot);
            r1.PEa1Tv1(timeStep, vTot);

            for (int iPass = 0; iPass < passes; iPass++) {
                double tStep = iPass == 0 ? tStep0 : (timeStep - tStep0);

                dr.Ev1Mv2(r1, r0);
                dr.ME(offset);
                boundary.nearestImage(dr);
                dv.Ev1Mv2(v1, v0);

                drTmp.E(dr);
                double rm0 = atom0.getType().rm();
                double rm1 = atom1.getType().rm();
                double rm2 = rm0 + rm1;
                omega = Math.sqrt(2 * w * rm2);
                cosomegat = Math.cos(omega * tStep);
                sinomegat = Math.sin(omega * tStep);
                drTmp.TE(cosomegat);
                drTmp.PEa1Tv1(sinomegat / omega, dv);
                drTmp.ME(dr);
                r0.PEa1Tv1(-rm0 / rm2, drTmp);
                r1.PEa1Tv1(+rm1 / rm2, drTmp);

                dvTmp.E(dv);
                dvTmp.TE(cosomegat);
                dvTmp.PEa1Tv1(-omega * sinomegat, dr);
                dvTmp.ME(dv);
                v0.PEa1Tv1(-rm0 / rm2, dvTmp);
                v1.PEa1Tv1(+rm1 / rm2, dvTmp);
                if (iPass == 0 && passes == 2) {
                    // randomize now
                    for (int i = 0; i < dv.getD(); i++) {
                        dvTmp.setX(i, random.nextGaussian());
                    }
                    dvTmp.TE(Math.sqrt(temperature * rm2));
                    v0.E(vTot);
                    v1.E(vTot);
                    v0.PEa1Tv1(-rm0 / rm2, dvTmp);
                    v1.PEa1Tv1(+rm1 / rm2, dvTmp);
                }
            }
        }

        eventManager.forcePrecomputed();

        forceSum.reset();
        //Compute forces on each atom
        potentialMaster.calculate(box, allAtoms, forceSum);

        eventManager.forceComputed();

        if (forceSum instanceof PotentialCalculationForcePressureSum) {
            pressureTensor.E(((PotentialCalculationForcePressureSum) forceSum).getPressureTensor());
        }

        //Finish integration step
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
//            System.out.println("force: "+((MyAgent)a.ia).force.toString());
            Vector velocity = a.getVelocity();
            workTensor.Ev1v2(velocity, velocity);
            workTensor.TE(((IAtom) a).getType().getMass());
            pressureTensor.PE(workTensor);
            if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                System.out.println("second " + a + " v=" + velocity + ", f=" + agentManager.getAgent(a).force);
            }
            velocity.PEa1Tv1(0.5 * timeStep * a.getType().rm(), agentManager.getAgent(a).force);  //p += f(new)*dt/2
        }

        pressureTensor.TE(1 / box.getBoundary().volume());

        if(isothermal) {
            doThermostatInternal();
        }
    }
}
