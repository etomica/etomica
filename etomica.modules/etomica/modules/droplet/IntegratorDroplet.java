package etomica.modules.droplet;

import java.io.Serializable;

import etomica.api.IAtomKinetic;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomList;
import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.api.IVectorMutable;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorMD;
import etomica.potential.PotentialCalculationForcePressureSum;
import etomica.potential.PotentialCalculationForceSum;
import etomica.space.ISpace;
import etomica.space.Tensor;
import etomica.util.Debug;

/**
 * Mesoscale integrator for Droplet module.
 * @author Andrew Schultz
 */
public class IntegratorDroplet extends IntegratorMD implements AgentSource {

    private static final long serialVersionUID = 2L;
    protected PotentialCalculationForceSum forceSum;;
    protected final IteratorDirective allAtoms;
    protected final Tensor pressureTensor;
    protected final Tensor workTensor, workTensor2;
    protected final Tensor identity;
    protected final IVectorMutable dr;

    protected AtomLeafAgentManager agentManager;

    public IntegratorDroplet(ISimulation sim, IPotentialMaster potentialMaster, ISpace _space) {
        this(potentialMaster, sim.getRandom(), 0.05, 1.0, _space);
    }

    public IntegratorDroplet(IPotentialMaster potentialMaster, IRandom random,
            double timeStep, double temperature, ISpace _space) {
        super(potentialMaster,random,timeStep,temperature, _space);
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

    public void setForceSum(PotentialCalculationForceSum pc){
        forceSum = pc;
        if(box != null){
            forceSum.setAgentManager(agentManager);
        }
        
    }

    public void setBox(IBox p) {
        if (box != null) {
            // allow agentManager to de-register itself as a BoxListener
            agentManager.dispose();
        }
        super.setBox(p);
        agentManager = new AtomLeafAgentManager(this,p);
        forceSum.setAgentManager(agentManager);
    }

//--------------------------------------------------------------
// steps all particles across time interval tStep

    // assumes one box
    public void doStepInternal() {
        super.doStepInternal();
        if (Debug.ON && Debug.DEBUG_NOW) {
            IAtomList pair = Debug.getAtoms(box);
            if (pair != null) {
                IVectorMutable dr = space.makeVector();
                dr.Ev1Mv2(((IAtomPositioned)pair.getAtom(1)).getPosition(), ((IAtomPositioned)pair.getAtom(0)).getPosition());
                System.out.println(pair+" dr "+dr);
            }
        }
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            MyAgent agent = (MyAgent)agentManager.getAgent((IAtomLeaf)a);
            agent.r0.E(((IAtomPositioned)a).getPosition());
            agent.rp.E(((IAtomPositioned)a).getPosition());
        }
        
        foo();

        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            MyAgent agent = (MyAgent)agentManager.getAgent((IAtomLeaf)a);
            IVectorMutable r = ((IAtomPositioned)a).getPosition();
            r.E(agent.r0);
            r.PEa1Tv1(0.5*timeStep, a.getVelocity());
            agent.rp.PEa1Tv1(timeStep/6.0, a.getVelocity());
//            if (((IAtomPositioned)a).getPosition().isNaN() || agent.rp.isNaN() || a.getVelocity().isNaN()) {
//                System.out.println("yoyo "+a+" "+ ((IAtomPositioned)a).getPosition());
//                System.out.println("rp "+agent.rp);
//                System.out.println("v "+a.getVelocity());
//                throw new RuntimeException("oops");
//            }
        }
        
        foo();
        
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            MyAgent agent = (MyAgent)agentManager.getAgent((IAtomLeaf)a);
            IVectorMutable r = ((IAtomPositioned)a).getPosition();
            r.E(agent.r0);
            r.PEa1Tv1(0.5*timeStep, a.getVelocity());
            agent.rp.PEa1Tv1(timeStep/3.0, a.getVelocity());
//            if (((IAtomPositioned)a).getPosition().isNaN() || agent.rp.isNaN() || a.getVelocity().isNaN()) {
//                System.out.println("bar "+a+" "+ ((IAtomPositioned)a).getPosition());
//                System.out.println("rp "+agent.rp);
//                System.out.println("v "+a.getVelocity());
//                throw new RuntimeException("oops");
//            }
        }
        
        foo();

        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            MyAgent agent = (MyAgent)agentManager.getAgent((IAtomLeaf)a);
            IVectorMutable r = ((IAtomPositioned)a).getPosition();
            r.E(agent.r0);
            r.PEa1Tv1(timeStep, a.getVelocity());
            agent.rp.PEa1Tv1(timeStep/3.0, a.getVelocity());
//            if (((IAtomPositioned)a).getPosition().isNaN() || agent.rp.isNaN() || a.getVelocity().isNaN()) {
//                System.out.println("foo "+a+" "+ ((IAtomPositioned)a).getPosition());
//                System.out.println("rp "+agent.rp);
//                System.out.println("v "+a.getVelocity());
//                throw new RuntimeException("oops");
//            }
        }
        
        foo();

        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            MyAgent agent = (MyAgent)agentManager.getAgent((IAtomLeaf)a);
            agent.rp.PEa1Tv1(timeStep/6.0, a.getVelocity());
            ((IAtomPositioned)a).getPosition().E(agent.rp);
//            if (((IAtomPositioned)a).getPosition().isNaN()) {
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
            IVectorMutable v = a.getVelocity();
            v.E(0);
        }
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            IVectorMutable v = a.getVelocity();
            dr.E(0);
            stokeslet(sp);
            MyAgent iAgent = (MyAgent)agentManager.getAgent((IAtomLeaf)a);
            dr.Ea1Tv1(dv,iAgent.force);
            workTensor.transform(dr);
            v.PE(dr);
            for (int jLeaf=iLeaf+1; jLeaf<nLeaf; jLeaf++) {
                IAtomKinetic aj = (IAtomKinetic)leafList.getAtom(jLeaf);
                MyAgent jAgent = (MyAgent)agentManager.getAgent((IAtomLeaf)aj);
                dr.Ev1Mv2(a.getPosition(),aj.getPosition());
                stokeslet(sp);
                // reuse as tensor * f
                dr.Ea1Tv1(dv,jAgent.force);
                workTensor.transform(dr);
                v.PE(dr);

                dr.Ea1Tv1(dv,iAgent.force);
                workTensor.transform(dr);
                IVectorMutable vj = aj.getVelocity();
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
    
    public void reset() throws ConfigurationOverlapException{
        if(!initialized) return;
        
        super.reset();
        if (Debug.ON && Debug.DEBUG_NOW) {
            IAtomList pair = Debug.getAtoms(box);
            if (pair != null) {
                IVectorMutable dr = space.makeVector();
                dr.Ev1Mv2(((IAtomPositioned)pair.getAtom(1)).getPosition(), ((IAtomPositioned)pair.getAtom(0)).getPosition());
                System.out.println(pair+" dr "+dr);
            }
        }

        forceSum.reset();
        potentialMaster.calculate(box, allAtoms, forceSum);
    }

//--------------------------------------------------------------
    
    public Class getAgentClass() {
        return MyAgent.class;
    }

    public final Object makeAgent(IAtomLeaf a) {
        return new MyAgent(space);
    }
    
    public void releaseAgent(Object agent, IAtomLeaf atom) {}
            
    public final static class MyAgent implements IntegratorBox.Forcible, Serializable {  //need public so to use with instanceof
        private static final long serialVersionUID = 1L;
        public IVectorMutable force;
        public IVectorMutable r0; // position at the beginning of the timestep
        public IVectorMutable rp;

        public MyAgent(ISpace space) {
            force = space.makeVector();
            r0 = space.makeVector();
            rp = space.makeVector();
        }
        
        public IVectorMutable force() {return force;}
    }
    
}
