// includes a main method

package etomica.integrator;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.exception.ConfigurationOverlapException;
import etomica.phase.Phase;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.ICoordinateKinetic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.systems.LJ;

/**
 * Gear 4th-order predictor-corrector integrator.
 *
 * @author Ed Maginn
 * @author David Kofke
 */
public class IntegratorGear4 extends IntegratorMD implements EtomicaElement, AgentSource {

    private final PotentialCalculationForceSum forceSum;
    private final IteratorDirective allAtoms = new IteratorDirective();
    protected final Space space;
    final Vector work1, work2;
    double zeta = 0.0;
    double chi = 0.0;
    double p1, p2, p3, p4;
    double c0, c2, c3, c4;
    
    static final double GEAR0 = 251./720.;
    static final double GEAR2 = 11./12.;
    static final double GEAR3 = 1./3.;
    static final double GEAR4 = 1./24.;

    protected Agent[] agents;
    protected AtomAgentManager agentManager;

    public IntegratorGear4(Simulation sim) {
        this(sim.potentialMaster,sim.space,sim.getDefaults().timeStep,sim.getDefaults().temperature);
    }
    
    public IntegratorGear4(PotentialMaster potentialMaster, Space space, 
            double timeStep, double temperature) {
        super(potentialMaster,timeStep,temperature);
        this.space = space;
        forceSum = new PotentialCalculationForceSum(space);
        work1 = space.makeVector();
        work2 = space.makeVector();
        //XXX this is totally wrong!  This should be based on the actual temperature and
        //potentials (steepness and depth) used.
        setTimeStep(new LJ().time().toSim(2.0));
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Molecular dynamics using Gear 4th-order predictor/corrector algorithm");
        return info;
    }

    public void setPhase(Phase p) {
        if (phase != null) {
            // allow agentManager to de-register itself as a PhaseListener
            agentManager.setPhase(null);
        }
        super.setPhase(p);
        agentManager = new AtomAgentManager(this,p);
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
        
    
    public void doStep() {
        
        predictor();
        calculateForces();
        corrector();
        
    }//end of doStep
    
    protected void calculateForces() {
        
        //Compute all forces
        atomIterator.reset();
        while(atomIterator.hasNext()) {   //zero forces on all atoms
            agents[atomIterator.nextAtom().getGlobalIndex()].force.E(0.0);
        }
        //Compute forces on each atom
        potential.calculate(phase, allAtoms, forceSum);
        
    }//end of calculateForces
    
    protected void corrector() {
        
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            AtomLeaf a = (AtomLeaf)atomIterator.nextAtom();
            Agent agent = agents[a.getGlobalIndex()];
            Vector r = a.coord.position();
            Vector v = ((ICoordinateKinetic)a.coord).velocity();
            work1.E(v);
            work1.PEa1Tv1(chi,r);
            work2.E(work1);
            work2.ME(agent.dr1);
            r.PEa1Tv1(c0, work2);
            agent.dr1.E(work1);
            agent.dr2.PEa1Tv1(c2,work2);
            agent.dr3.PEa1Tv1(c3,work2);
            agent.dr4.PEa1Tv1(c4,work2);
            
            work1.Ea1Tv1(((AtomTypeLeaf)a.type).rm(),agent.force);
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
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            AtomLeaf a = (AtomLeaf)atomIterator.nextAtom();
            Agent agent = agents[a.getGlobalIndex()];
            Vector r = a.coord.position();
            Vector v = ((ICoordinateKinetic)a.coord).velocity();
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
        agents = (Agent[])agentManager.getAgents();
        calculateForces();
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            AtomLeaf a = (AtomLeaf)atomIterator.nextAtom();
            Agent agent = agents[a.getGlobalIndex()];
            agent.dr1.E(((ICoordinateKinetic)a.coord).velocity());
            agent.dr2.Ea1Tv1(((AtomTypeLeaf)a.type).rm(),agent.force);
            agent.dr3.E(0.0);
            agent.dr4.E(0.0);
            agent.dv1.Ea1Tv1(((AtomTypeLeaf)a.type).rm(),agent.force);
            agent.dv2.E(0.0);
            agent.dv3.E(0.0);
            agent.dv4.E(0.0);
        }
        super.reset();
    }

    public Class getAgentClass() {
        return Agent.class;
    }
    
    public Object makeAgent(Atom a) {
        return new Agent(space,a);
    }
    
    public void releaseAgent(Object agent, Atom atom) {}
            
    public static class Agent implements IntegratorPhase.Forcible {  //need public so to use with instanceof
        public Atom atom;
        public Vector force;
        public Vector dr1, dr2, dr3, dr4;
        public Vector dv1, dv2, dv3, dv4;

        public Agent(Space space, Atom a) {
            atom = a;
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
        
        public Vector force() {return force;}
    }
}
