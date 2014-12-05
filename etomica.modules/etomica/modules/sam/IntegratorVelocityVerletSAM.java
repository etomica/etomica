/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.sam;

import etomica.api.IAtom;
import etomica.api.IAtomKinetic;
import etomica.api.IAtomList;
import etomica.api.IAtomType;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.IVectorMutable;
import etomica.atom.AtomSetSinglet;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.PotentialCalculationForcePressureSum;
import etomica.space.ISpace;
import etomica.util.Debug;

/**
 * Special integrator subclass that has the ability to prevent sulfur atoms
 * from moving in the Y direction.
 *
 * @author Andrew Schultz
 */
public class IntegratorVelocityVerletSAM extends IntegratorVelocityVerlet {

    public IntegratorVelocityVerletSAM(IPotentialMaster potentialMaster,
            IRandom random, double timeStep, double temperature, ISpace _space) {
        super(potentialMaster, random, timeStep, temperature, _space);
    }
    
    public void setSulfurType(IAtomType newSulfurType) {
        sulfurType = newSulfurType;
    }
    
    public void doStepInternal() {
        super.doStepInternal();
        if (Debug.ON && Debug.DEBUG_NOW) {
            IAtomList pair = Debug.getAtoms(box);
            if (pair != null) {
                IVectorMutable dr = space.makeVector();
                dr.Ev1Mv2(pair.getAtom(1).getPosition(), pair.getAtom(0).getPosition());
                System.out.println(pair+" dr "+dr);
            }
        }
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            MyAgent agent = (MyAgent)agentManager.getAgent((IAtom)a);
            IVectorMutable r = a.getPosition();
            IVectorMutable v = a.getVelocity();
            if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                System.out.println("first "+a+" r="+r+", v="+v+", f="+agent.force);
            }
            v.PEa1Tv1(0.5*timeStep*((IAtom)a).getType().rm(),agent.force);  // p += f(old)*dt/2
            if (((IAtom)a).getType() == sulfurType) {
                // sulfur isn't allowed to move in the Y direction
                v.setX(1, 0);
            }
            r.PEa1Tv1(timeStep,v);         // r += p*dt/m
        }

        forceSum.reset();
        //Compute forces on each atom
        potentialMaster.calculate(box, allAtoms, forceSum);
        
        if(forceSum instanceof PotentialCalculationForcePressureSum){
            pressureTensor.E(((PotentialCalculationForcePressureSum)forceSum).getPressureTensor());
        }
        
        //Finish integration step
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
//            System.out.println("force: "+((MyAgent)a.ia).force.toString());
            IVectorMutable velocity = a.getVelocity();
            workTensor.Ev1v2(velocity,velocity);
            workTensor.TE(((IAtom)a).getType().getMass());
            pressureTensor.PE(workTensor);
            if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                System.out.println("second "+a+" v="+velocity+", f="+((MyAgent)agentManager.getAgent(a)).force);
            }
            velocity.PEa1Tv1(0.5*timeStep*((IAtom)a).getType().rm(),((MyAgent)agentManager.getAgent(a)).force);  //p += f(new)*dt/2
            if (((IAtom)a).getType() == sulfurType) {
                velocity.setX(1, 0);
            }
        }
        
        pressureTensor.TE(1/box.getBoundary().volume());

        if(isothermal) {
            doThermostatInternal();
        }
    }
    
    protected IAtomType sulfurType;
}
