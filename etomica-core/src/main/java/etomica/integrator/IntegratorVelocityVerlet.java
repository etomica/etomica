/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.atom.*;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.box.Box;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForcePressureSum;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.util.Debug;
import etomica.util.random.IRandom;

public class IntegratorVelocityVerlet extends IntegratorMD implements AgentSource<Vector> {

    protected PotentialCalculationForceSum forceSum;;
    protected final IteratorDirective allAtoms;
    protected final Tensor pressureTensor;
    protected final Tensor workTensor;

    protected AtomLeafAgentManager<Vector> agentManager;

    public IntegratorVelocityVerlet(Simulation sim, PotentialMaster potentialMaster, Box box) {
        this(potentialMaster, sim.getRandom(), 0.05, 1.0, box);
    }

    public IntegratorVelocityVerlet(PotentialMaster potentialMaster, IRandom random,
                                    double timeStep, double temperature, Box box) {
        super(potentialMaster,random,timeStep,temperature, box);
        // if you're motivated to throw away information earlier, you can use 
        // PotentialCalculationForceSum instead.
        forceSum = new PotentialCalculationForcePressureSum(space);
        allAtoms = new IteratorDirective();
        // allAtoms is used only for the force calculation, which has no LRC
        // but we're also calculating the pressure tensor, which does have LRC.
        // things deal with this OK.
        allAtoms.setIncludeLrc(true);
        pressureTensor = space.makeTensor();
        workTensor = space.makeTensor();
        agentManager = new AtomLeafAgentManager<>(this, box);
        forceSum.setAgentManager(agentManager);
    }

    public PotentialCalculationForceSum getForceSum() {
        return forceSum;
    }

    public AtomLeafAgentManager<Vector> getAgentManager() {
        return agentManager;
    }

    public void setForceSum(PotentialCalculationForceSum pc){
        forceSum = pc;
        if(box != null){
            forceSum.setAgentManager(agentManager);
        }
        
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
                System.out.println(pair+" dr "+dr);
            }
        }
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.size();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.get(iLeaf);
            Vector force = agentManager.getAgent(a);
            Vector r = a.getPosition();
            Vector v = a.getVelocity();
            if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                System.out.println("first "+a+" r="+r+", v="+v+", f="+force);
            }
            v.PEa1Tv1(0.5*timeStep* a.getType().rm(),force);  // p += f(old)*dt/2
            r.PEa1Tv1(timeStep,v);         // r += p*dt/m
        }
    
        eventManager.forcePrecomputed();
    
        forceSum.reset();
        //Compute forces on each atom
        potentialMaster.calculate(box, allAtoms, forceSum);

        eventManager.forceComputed();

        if(forceSum instanceof PotentialCalculationForcePressureSum){
            pressureTensor.E(((PotentialCalculationForcePressureSum)forceSum).getPressureTensor());
        }

        //Finish integration step
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.get(iLeaf);
//            System.out.println("force: "+((MyAgent)a.ia).force.toString());
            Vector velocity = a.getVelocity();
            workTensor.Ev1v2(velocity,velocity);
            workTensor.TE(a.getType().getMass());
            pressureTensor.PE(workTensor);
            if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                System.out.println("second "+a+" v="+velocity+", f="+ agentManager.getAgent(a));
            }
            velocity.PEa1Tv1(0.5*timeStep* a.getType().rm(), agentManager.getAgent(a));  //p += f(new)*dt/2
        }

        pressureTensor.TE(1/box.getBoundary().volume());

        if(isothermal) {
            doThermostatInternal();
        }
    }

    /**
     * Returns the pressure tensor based on the forces calculated during the
     * last time step.
     */
    public Tensor getPressureTensor() {
        return pressureTensor;
    }
    
    public void reset() {
        super.reset();

        if (Debug.ON && Debug.DEBUG_NOW) {
            IAtomList pair = Debug.getAtoms(box);
            if (pair != null) {
                Vector dr = space.makeVector();
                dr.Ev1Mv2(pair.get(1).getPosition(), pair.get(0).getPosition());
                System.out.println(pair+" dr "+dr);
            }
        }

        precomputeForce();
    }

    public void precomputeForce() {
        eventManager.forcePrecomputed();

        forceSum.reset();
        potentialMaster.calculate(box, allAtoms, forceSum);

        eventManager.forceComputed();
    }

//--------------------------------------------------------------

    public Vector makeAgent(IAtom a, Box agentBox) {
        return space.makeVector();
    }

    public void releaseAgent(Vector agent, IAtom atom, Box agentBox) {
    }

    public void postRestore() {
        super.postRestore();
        precomputeForce();
    }
}
