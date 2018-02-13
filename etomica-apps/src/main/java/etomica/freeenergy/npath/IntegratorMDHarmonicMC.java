/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.freeenergy.npath;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.PotentialMaster;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Created by andrew on 4/30/17.
 */
public class IntegratorMDHarmonicMC extends IntegratorVelocityVerlet {

    protected final Vector dr, drTmp, dv, dvTmp;
    protected final Vector fTot, vTot;
    protected P1ImageHarmonic p1;
    protected Boundary boundary;
    protected Vector[] drAll;
    protected boolean firstTrial = true;
    protected int nAttempted;
    protected double chiSum;

    public IntegratorMDHarmonicMC(PotentialMaster potentialMaster, IRandom random, double timeStep, double temperature, Space space, Box box) {
        super(potentialMaster, random, timeStep, temperature, space, box);
        dr = space.makeVector();
        drTmp = space.makeVector();
        dv = space.makeVector();
        dvTmp = space.makeVector();
        fTot = space.makeVector();
        vTot = space.makeVector();
    }

    public void setP1Harmonic(P1ImageHarmonic p1) {
        this.p1 = p1;
    }

    public void setBox(Box box) {
        super.setBox(box);
        boundary = box.getBoundary();
        if (drAll == null || drAll.length != box.getLeafList().getAtomCount()) {
            drAll = space.makeVectorArray(box.getLeafList().getAtomCount());
            firstTrial = true;
        }
    }

    protected void doStepInternal() {
        // from IntegratorMD
        currentTime += timeStep;
        // IntegratorVelocityVerlet code
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            int iLeaf1 = p1.getPartner(iLeaf);
            if (iLeaf1 < iLeaf) {
                continue;
            }
            IAtomKinetic a = (IAtomKinetic) leafList.getAtom(iLeaf);
            Vector force = agentManager.getAgent(a);
            Vector v = a.getVelocity();
            IAtomKinetic a1 = (IAtomKinetic) leafList.getAtom(iLeaf1);
            Vector force1 = agentManager.getAgent(a1);
            Vector v1 = a1.getVelocity();
            fTot.Ev1Pv2(force, force1);
            dv.Ea1Tv1(0.25 * timeStep * a.getType().rm(), fTot);
            v.PE(dv);
            v1.PE(dv);

            vTot.Ev1Pv2(v1, v);
            vTot.TE(0.5);
            a.getPosition().PEa1Tv1(timeStep, vTot);
            a1.getPosition().PEa1Tv1(timeStep, vTot);
        }

        Vector offset = p1.getOffset();
        double w = p1.getW();
        // p = exp(-wx^2) = exp(-x^2/2sigma^2)
        // w = 1/2sigma^2
        double sigma = Math.sqrt(0.5 * temperature / w);
        double u0 = meterPE.getDataAsScalar();
        double du = 0;
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            int iLeaf1 = p1.getPartner(iLeaf);
            if (iLeaf1 < iLeaf) {
                continue;
            }
            IAtom atom0 = leafList.getAtom(iLeaf);
            IAtom atom1 = leafList.getAtom(iLeaf1);
            Vector r0 = atom0.getPosition();
            Vector r1 = atom1.getPosition();

            dr.Ev1Mv2(r1, r0);
            dr.ME(offset);
            boundary.nearestImage(dr);
            drAll[iLeaf].E(dr);
            du -= w * dr.squared();

            // f = -2 w (x-x0) + f_iron
            // f = -2 w ((x-x0) - f_iron/2w)
            // f = -2 w (x - (x0 + f_iron/2w))

            // (1/2sigma^2) = w/T
            // sigma = sqrt(2T/w)
            for (int k = 0; k < dr.getD(); k++) {
                drTmp.setX(k, random.nextGaussian() * sigma);
            }
            du += w * drTmp.squared();
//            drTmp.PEa1Tv1(0.5/w, df);
            drTmp.ME(dr);
            r0.PEa1Tv1(-0.5, drTmp);
            r1.PEa1Tv1(+0.5, drTmp);
        }
        double u = meterPE.getDataAsScalar();
        double x = u - u0 - du;
        double chi = x < 0 ? 1 : Math.exp(-x / temperature);
        // we always accept the first move because it's really difficult
        if (!firstTrial && chi < 1 && chi < random.nextDouble()) {
            for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
                int iLeaf1 = p1.getPartner(iLeaf);
                if (iLeaf1 < iLeaf) {
                    continue;
                }
                IAtom atom0 = leafList.getAtom(iLeaf);
                IAtom atom1 = leafList.getAtom(iLeaf1);
                Vector r0 = atom0.getPosition();
                Vector r1 = atom1.getPosition();

                dr.Ev1Mv2(r1, r0);
                dr.ME(offset);
                boundary.nearestImage(dr);
                drAll[iLeaf].ME(dr);
                r0.PEa1Tv1(-0.5, drAll[iLeaf]);
                r1.PEa1Tv1(+0.5, drAll[iLeaf]);
            }
        }
        chiSum += chi;
        nAttempted++;
        firstTrial = false;

        eventManager.forcePrecomputed();

        forceSum.reset();
        //Compute forces on each atom
        potentialMaster.calculate(box, allAtoms, forceSum);

        eventManager.forceComputed();

        //Finish integration step
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            int iLeaf1 = p1.getPartner(iLeaf);
            if (iLeaf1 < iLeaf) {
                continue;
            }
            IAtomKinetic a = (IAtomKinetic) leafList.getAtom(iLeaf);
            IAtomKinetic a1 = (IAtomKinetic) leafList.getAtom(iLeaf1);
            fTot.Ev1Pv2(agentManager.getAgent(a), agentManager.getAgent(a1));
            dv.Ea1Tv1(0.25 * timeStep * a.getType().rm(), fTot);
            a.getVelocity().PE(dv);
            a1.getVelocity().PE(dv);
        }

        if (isothermal) {
            doThermostatInternal();
        }
    }

    public void resetAcceptance() {
        chiSum = 0;
        nAttempted = 0;
    }

    public double getAcceptanceProbability() {
        return chiSum / nAttempted;
    }
}
