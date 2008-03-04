// includes a main method

package etomica.integrator;

import etomica.EtomicaInfo;
import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.api.IVector;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.IAtomKinetic;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.exception.ConfigurationOverlapException;
import etomica.potential.PotentialCalculationForceSum;
import etomica.space.Space;

/**
 * Gear 4th-order predictor-corrector integrator.
 *
 * @author Ed Maginn
 * @author David Kofke
 */
public class IntegratorGear4 extends IntegratorMD implements AgentSource {

    private static final long serialVersionUID = 1L;
    private final PotentialCalculationForceSum forceSum;
    private final IteratorDirective allAtoms;
    final IVector work1, work2;
    double zeta = 0.0;
    double chi = 0.0;
    double p1, p2, p3, p4;
    double c0, c2, c3, c4;
    
    static final double GEAR0 = 251./720.;
    static final double GEAR2 = 11./12.;
    static final double GEAR3 = 1./3.;
    static final double GEAR4 = 1./24.;

    protected AtomLeafAgentManager agentManager;

    public IntegratorGear4(ISimulation sim, IPotentialMaster potentialMaster, Space _space) {
        this(potentialMaster, sim.getRandom(), 0.05, 1.0, _space);
    }
    
    public IntegratorGear4(IPotentialMaster potentialMaster, IRandom random, 
            double timeStep, double temperature, Space _space) {
        super(potentialMaster,random,timeStep,temperature, _space);
        forceSum = new PotentialCalculationForceSum();
        allAtoms = new IteratorDirective();
        // allAtoms is used only for the force calculation, which has no LRC
        allAtoms.setIncludeLrc(false);
        work1 = potentialMaster.getSpace().makeVector();
        work2 = potentialMaster.getSpace().makeVector();

        setTimeStep(timeStep);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Molecular dynamics using Gear 4th-order predictor/corrector algorithm");
        return info;
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

    public void setTimeStep(double dt) {
        super.setTimeStep(dt);
        p1 = dt;
        p2 = p1 * dt / 2.0;
        p3 = p2 * dt / 3.0;
        p4 = p3 * dt / 4.0;
        c0 = GEAR0 * p1;
        c2 = GEAR2 * p1 / p2;
        c3 = GEAR3 * p1 / p3;
        c4 = GEAR4 * p1 / p4;
    }
        
    
    public void doStepInternal() {
        super.doStepInternal();
        
        predictor();
        calculateForces();
        corrector();
        
    }//end of doStep
    
    protected void calculateForces() {
        //Compute all forces

        //zero forces on all atoms
        forceSum.reset();
        //Compute forces on each atom
        potential.calculate(box, allAtoms, forceSum);
        
    }//end of calculateForces
    
    protected void corrector() {
        
        IAtomSet leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            Agent agent = (Agent)agentManager.getAgent(a);
            IVector r = a.getPosition();
            IVector v = a.getVelocity();
            work1.E(v);
            work1.PEa1Tv1(chi,r);
            work2.E(work1);
            work2.ME(agent.dr1);
            r.PEa1Tv1(c0, work2);
            if (r.isNaN()) {
                throw new RuntimeException("oops "+a+" "+r+" "+work1+" "+work2+" "+chi+" "+c0);
            }
            agent.dr1.E(work1);
            agent.dr2.PEa1Tv1(c2,work2);
            agent.dr3.PEa1Tv1(c3,work2);
            agent.dr4.PEa1Tv1(c4,work2);
            
            work1.Ea1Tv1(((AtomTypeLeaf)a.getType()).rm(),agent.force);
            work1.PEa1Tv1(-(zeta+chi),v);
            work2.E(work1);
            work2.ME(agent.dv1);
            v.PEa1Tv1(c0,work2);
            agent.dv1.E(work1);
            agent.dv2.PEa1Tv1(c2,work2);
            agent.dv3.PEa1Tv1(c3,work2);
            agent.dv4.PEa1Tv1(c4,work2);
        }
    }//end of corrector
        
    protected void predictor() {
        IAtomSet leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            Agent agent = (Agent)agentManager.getAgent(a);
            IVector r = a.getPosition();
            IVector v = a.getVelocity();
            r.PEa1Tv1(p1, agent.dr1);
            r.PEa1Tv1(p2, agent.dr2);
            r.PEa1Tv1(p3, agent.dr3);
            r.PEa1Tv1(p4, agent.dr4);
            
            agent.dr1.PEa1Tv1(p1, agent.dr2);
            agent.dr1.PEa1Tv1(p2, agent.dr3);
            agent.dr1.PEa1Tv1(p3, agent.dr4);
            
            agent.dr2.PEa1Tv1(p1, agent.dr3);
            agent.dr2.PEa1Tv1(p2, agent.dr4);
            
            agent.dr3.PEa1Tv1(p1, agent.dr4);
            
            v.PEa1Tv1(p1, agent.dv1);
            v.PEa1Tv1(p2, agent.dv2);
            v.PEa1Tv1(p3, agent.dv3);
            v.PEa1Tv1(p4, agent.dv4);
            
            agent.dv1.PEa1Tv1(p1, agent.dv2);
            agent.dv1.PEa1Tv1(p2, agent.dv3);
            agent.dv1.PEa1Tv1(p3, agent.dv4);
            
            agent.dv2.PEa1Tv1(p1, agent.dv3);
            agent.dv2.PEa1Tv1(p2, agent.dv4);
            
            agent.dv3.PEa1Tv1(p1, agent.dv4);
        }
    }

    public void reset() throws ConfigurationOverlapException {
        calculateForces();
        IAtomSet leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            Agent agent = (Agent)agentManager.getAgent(a);
            agent.dr1.E(a.getVelocity());
            agent.dr2.Ea1Tv1(((AtomTypeLeaf)a.getType()).rm(),agent.force);
            agent.dr3.E(0.0);
            agent.dr4.E(0.0);
            agent.dv1.Ea1Tv1(((AtomTypeLeaf)a.getType()).rm(),agent.force);
            agent.dv2.E(0.0);
            agent.dv3.E(0.0);
            agent.dv4.E(0.0);
        }
        super.reset();
    }

    public Class getAgentClass() {
        return Agent.class;
    }
    
    public Object makeAgent(IAtom a) {
        return new Agent(space);
    }
    
    public void releaseAgent(Object agent, IAtom atom) {}
            
    public static class Agent implements IntegratorBox.Forcible {  //need public so to use with instanceof
        public IVector force;
        public IVector dr1, dr2, dr3, dr4;
        public IVector dv1, dv2, dv3, dv4;

        public Agent(Space space) {
            force = space.makeVector();
            dr1 = space.makeVector();
            dr2 = space.makeVector();
            dr3 = space.makeVector();
            dr4 = space.makeVector();
            dv1 = space.makeVector();
            dv2 = space.makeVector();
            dv3 = space.makeVector();
            dv4 = space.makeVector();
        }
        
        public IVector force() {return force;}
    }
}
