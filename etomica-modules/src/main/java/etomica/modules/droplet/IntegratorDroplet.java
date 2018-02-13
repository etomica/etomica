/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.droplet;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorMD;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForcePressureSum;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.util.Debug;
import etomica.util.random.IRandom;

import java.io.Serializable;

/**
 * Mesoscale integrator for Droplet module.
 * @author Andrew Schultz
 */
public class IntegratorDroplet extends IntegratorMD implements AgentSource<IntegratorDroplet.MyAgent> {

    private static final long serialVersionUID = 2L;
    protected PotentialCalculationForceSum forceSum;;
    protected final IteratorDirective allAtoms;
    protected final Tensor pressureTensor;
    protected final Tensor workTensor, workTensor2;
    protected final Tensor identity;
    protected final Vector dr;

    protected AtomLeafAgentManager<MyAgent> agentManager;
    private AtomLeafAgentManager<Vector> forces;

    public IntegratorDroplet(Simulation sim, PotentialMaster potentialMaster, Space _space, Box box) {
        this(potentialMaster, sim.getRandom(), 0.05, 1.0, _space, box);
    }

    public IntegratorDroplet(PotentialMaster potentialMaster, IRandom random,
                             double timeStep, double temperature, Space _space, Box box) {
        super(potentialMaster,random,timeStep,temperature, _space, box);
        // if you're motivated to throw away information earlier, you can use 
        // PotentialCalculationForceSum instead.
        forceSum = new PotentialCalculationForcePressureSum(space);
        allAtoms = new IteratorDirective();
        // allAtoms is used only for the force calculation, which has no LRC
        // but we're also calculating the pressure tensor, which does have LRC.
        // things deal with this OK.
        allAtoms.setIncludeLrc(true);
        pressureTensor = space.makeTensor();
        identity = space.makeTensor();
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                identity.setComponent(i,j,i==j?1:0);
            }
        }
        workTensor = space.makeTensor();
        workTensor2 = space.makeTensor();
        dr = space.makeVector();
    }
    
    public void setBox(Box box) {
        if (this.box != null) {
            // allow agentManager to de-register itself as a BoxListener
            agentManager.dispose();
        }
        super.setBox(box);
        agentManager = new AtomLeafAgentManager<MyAgent>(this, box);
        forces = new AtomLeafAgentManager<>(a -> space.makeVector(), box);
        forceSum.setAgentManager(forces);
    }

//--------------------------------------------------------------
// steps all particles across time interval tStep

    // assumes one box
    protected void doStepInternal() {
        super.doStepInternal();

        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            MyAgent agent = agentManager.getAgent(a);
            agent.r0.E(a.getPosition());
            agent.rp.E(a.getPosition());
        }
        
        foo();

        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            MyAgent agent = agentManager.getAgent(a);
            Vector r = a.getPosition();
            r.E(agent.r0);
            r.PEa1Tv1(0.5*timeStep, a.getVelocity());
            agent.rp.PEa1Tv1(timeStep/6.0, a.getVelocity());
//            if (a).getPosition().isNaN() || agent.rp.isNaN() || a.getVelocity().isNaN()) {
//                System.out.println("yoyo "+a+" "+ a).getPosition());
//                System.out.println("rp "+agent.rp);
//                System.out.println("v "+a.getVelocity());
//                throw new RuntimeException("oops");
//            }
        }
        
        foo();
        
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            MyAgent agent = agentManager.getAgent(a);
            Vector r = a.getPosition();
            r.E(agent.r0);
            r.PEa1Tv1(0.5*timeStep, a.getVelocity());
            agent.rp.PEa1Tv1(timeStep/3.0, a.getVelocity());
//            if (a).getPosition().isNaN() || agent.rp.isNaN() || a.getVelocity().isNaN()) {
//                System.out.println("bar "+a+" "+ a).getPosition());
//                System.out.println("rp "+agent.rp);
//                System.out.println("v "+a.getVelocity());
//                throw new RuntimeException("oops");
//            }
        }
        
        foo();

        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            MyAgent agent = agentManager.getAgent(a);
            Vector r = a.getPosition();
            r.E(agent.r0);
            r.PEa1Tv1(timeStep, a.getVelocity());
            agent.rp.PEa1Tv1(timeStep/3.0, a.getVelocity());
//            if (a).getPosition().isNaN() || agent.rp.isNaN() || a.getVelocity().isNaN()) {
//                System.out.println("foo "+a+" "+ a).getPosition());
//                System.out.println("rp "+agent.rp);
//                System.out.println("v "+a.getVelocity());
//                throw new RuntimeException("oops");
//            }
        }
        
        foo();

        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            MyAgent agent = agentManager.getAgent(a);
            agent.rp.PEa1Tv1(timeStep/6.0, a.getVelocity());
            a.getPosition().E(agent.rp);
//            if (a).getPosition().isNaN()) {
//                System.out.println("rp "+agent.rp);
//                System.out.println("v "+a.getVelocity());
//                throw new RuntimeException("oops");
//            }
//            double vol = 4.0/3.0*Math.PI;
//            double dv = vol / box.getMoleculeList().getMoleculeCount();
//            agent.force.TE(dv);
//            System.out.println(a+" "+agent.rp+" "+a.getVelocity()+" "+agent.force);
        }
//        System.out.println(System.currentTimeMillis()/1000);
    }
    
    protected void foo() {
        forceSum.reset();
        //Compute forces on each atom
        potentialMaster.calculate(box, allAtoms, forceSum);
        double vol = 4.0/3.0*Math.PI;
        double dv = vol / box.getMoleculeList().getMoleculeCount();
        double sp = 2.0*Math.pow(dv,1.0/3.0);

        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            Vector v = a.getVelocity();
            v.E(0);
        }
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            Vector v = a.getVelocity();
            dr.E(0);
            stokeslet(sp);
            Vector iForce = forces.getAgent(a);
            dr.Ea1Tv1(dv,iForce);
            workTensor.transform(dr);
            v.PE(dr);
            for (int jLeaf=iLeaf+1; jLeaf<nLeaf; jLeaf++) {
                IAtomKinetic aj = (IAtomKinetic)leafList.getAtom(jLeaf);
                Vector jForce = forces.getAgent(aj);
                dr.Ev1Mv2(a.getPosition(),aj.getPosition());
                stokeslet(sp);
                // reuse as tensor * f
                dr.Ea1Tv1(dv,jForce);
                workTensor.transform(dr);
                v.PE(dr);

                dr.Ea1Tv1(dv,iForce);
                workTensor.transform(dr);
                Vector vj = aj.getVelocity();
                vj.PE(dr);
            }
        }
    }
    
    protected void stokeslet(double sp) {
        double r2 = dr.squared();
        double d = Math.sqrt(r2);
        double c1 = (0.5 - i(d/sp,2)) / (2 * Math.PI * sp);
        if (d < 1.e-9) {
            workTensor.E(identity);
            workTensor.TE(c1);
            return;
        }
        double c0 = 3*i(d/sp,3);
        double c2 = 0.5*i(d/sp,5) * sp*sp;
        workTensor2.Ev1v2(dr,dr);
        workTensor2.TE(1.0/(d*d));
        workTensor.E(workTensor2);
        workTensor.PE(identity);
        workTensor.TE(c0*0.125/(Math.PI*d));
        workTensor2.TE(-3.0);
        workTensor2.PE(identity);
        workTensor2.TE(c2*0.25/(Math.PI*d*d*d));
        workTensor.PE(workTensor2);
        workTensor2.E(identity);
        workTensor2.TE(c1);
        workTensor.PE(workTensor2);
    }
    
    protected double i(double x, int p) {
        if (x <= 0) {
            return 0;
        }
        if (x >= 1) {
            return 1.0/(p);
        }
        return Math.pow(x,p)/(p);
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
                dr.Ev1Mv2(pair.getAtom(1).getPosition(), pair.getAtom(0).getPosition());
                System.out.println(pair+" dr "+dr);
            }
        }

        forceSum.reset();
        potentialMaster.calculate(box, allAtoms, forceSum);
    }

//--------------------------------------------------------------

    public final MyAgent makeAgent(IAtom a, Box agentBox) {
        return new MyAgent(space);
    }
    
    public void releaseAgent(MyAgent agent, IAtom atom, Box agentBox) {}
            
    public final static class MyAgent {  //need public so to use with instanceof
        public Vector r0; // position at the beginning of the timestep
        public Vector rp;

        public MyAgent(Space space) {
            r0 = space.makeVector();
            rp = space.makeVector();
        }
    }
    
}
