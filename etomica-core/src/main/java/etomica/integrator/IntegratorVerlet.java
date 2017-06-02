/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.api.*;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.potential.PotentialCalculationForcePressureSum;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space.Tensor;

public final class IntegratorVerlet extends IntegratorMD implements AgentSource<IntegratorVerlet.Agent> {

    private static final long serialVersionUID = 1L;
    protected final PotentialCalculationForcePressureSum forceSum;
    private final IteratorDirective allAtoms;
    private double t2;
    protected final Tensor pressureTensor;
    protected final Tensor workTensor;

    Vector work;

    protected AtomLeafAgentManager<Agent> agentManager;

    public IntegratorVerlet(Simulation sim, PotentialMaster potentialMaster, Space _space) {
        this(potentialMaster, sim.getRandom(), 0.05, 1.0, _space);
    }
    
    public IntegratorVerlet(PotentialMaster potentialMaster, IRandom random,
                            double timeStep, double temperature, Space _space) {
        super(potentialMaster,random,timeStep,temperature, _space);
        // if you're motivated to throw away information earlier, you can use 
        // PotentialCalculationForceSum instead.
        forceSum = new PotentialCalculationForcePressureSum(space);
        allAtoms = new IteratorDirective();
        // but we're also calculating the pressure tensor, which does have LRC.
        // things deal with this OK.
        allAtoms.setIncludeLrc(true);
        work = space.makeVector();
        
        pressureTensor = space.makeTensor();
        workTensor = space.makeTensor();
    }

    public final void setTimeStep(double t) {
        super.setTimeStep(t);
        t2 = timeStep * timeStep;
    }
          
    public void setBox(Box p) {
        if (box != null) {
            // allow agentManager to de-register itself as a BoxListener
            agentManager.dispose();
        }
        super.setBox(p);
        agentManager = new AtomLeafAgentManager<Agent>(this,p,Agent.class);
        forceSum.setAgentManager(agentManager);
    }
    
//--------------------------------------------------------------
// steps all particles across time interval tStep

    public void doStepInternal() {
        super.doStepInternal();
        //Compute forces on each atom
        forceSum.reset();
        potentialMaster.calculate(box, allAtoms, forceSum);
        pressureTensor.E(forceSum.getPressureTensor());

        //take step
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            pressureTensor.E(forceSum.getPressureTensor());
            Vector v = a.getVelocity();
            workTensor.Ev1v2(v,v);
            workTensor.TE(a.getType().getMass());
            pressureTensor.PE(workTensor);
            
            Agent agent = agentManager.getAgent(a);
            Vector r = a.getPosition();
            work.E(r);
            r.PE(agent.rMrLast);
            agent.force.TE(a.getType().rm()*t2);
            r.PE(agent.force);
            agent.rMrLast.E(r);
            agent.rMrLast.ME(work);
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
        updateMrLast();
    }

    protected void updateMrLast() {
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            Agent agent = agentManager.getAgent(a);
            agent.rMrLast.Ea1Tv1(timeStep,a.getVelocity());//06/13/03 removed minus sign before timeStep
        }
    }

    /**
     * Updates MrLast appropriately after randomizing momenta
     * as part of the Andersen thermostat.
     */
    protected void randomizeMomenta() {
        super.randomizeMomenta();
        // super.randomizeMomenta changes the velocities, so we need to
        // recalculate hypothetical old positions
        updateMrLast();
    }

    /**
     * Updates MrLast appropriately after randomizing momentum
     * as part of the Andersen thermostat.
     */
    protected void randomizeMomentum(IAtomKinetic atom) {
        super.randomizeMomentum(atom);
        Agent agent = agentManager.getAgent(atom);
        agent.rMrLast.Ea1Tv1(timeStep,atom.getVelocity());//06/13/03 removed minus sign before timeStep
    }
    
    public final Agent makeAgent(IAtom a, Box agentBox) {
        return new Agent(space);
    }
    
    public void releaseAgent(Agent agent, IAtom atom, Box agentBox) {}
            
	public final static class Agent implements IntegratorBox.Forcible {  //need public so to use with instanceof
        public Vector force;
        public Vector rMrLast;  //r - rLast

        public Agent(Space space) {
            force = space.makeVector();
            rMrLast = space.makeVector();
        }
        
        public Vector force() {return force;}
    }//end of Agent
    
}//end of IntegratorVerlet

