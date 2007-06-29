package etomica.integrator;

import etomica.EtomicaInfo;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.exception.ConfigurationOverlapException;
import etomica.box.Box;
import etomica.potential.PotentialCalculationForcePressureSum;
import etomica.potential.PotentialMaster;
import etomica.simulation.ISimulation;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.util.IRandom;

public final class IntegratorVerlet extends IntegratorMD implements AgentSource {

    private static final long serialVersionUID = 1L;
    protected final PotentialCalculationForcePressureSum forceSum;
    private final IteratorDirective allAtoms;
    private double t2;
    protected final Tensor pressureTensor;
    protected final Tensor workTensor;

    IVector work;

    protected AtomLeafAgentManager agentManager;

    public IntegratorVerlet(ISimulation sim, PotentialMaster potentialMaster) {
        this(potentialMaster, sim.getRandom(), 0.05, 1.0);
    }
    
    public IntegratorVerlet(PotentialMaster potentialMaster, IRandom random, 
            double timeStep, double temperature) {
        super(potentialMaster,random,timeStep,temperature);
        // if you're motivated to throw away information earlier, you can use 
        // PotentialCalculationForceSum instead.
        forceSum = new PotentialCalculationForcePressureSum(potentialMaster.getSpace());
        allAtoms = new IteratorDirective();
        // but we're also calculating the pressure tensor, which does have LRC.
        // things deal with this OK.
        allAtoms.setIncludeLrc(true);
        work = potentialMaster.getSpace().makeVector();
        
        pressureTensor = potentialMaster.getSpace().makeTensor();
        workTensor = potentialMaster.getSpace().makeTensor();
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Molecular dynamics using basic Verlet algorithm");
        return info;
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
        agentManager = new AtomLeafAgentManager(this,p);
        forceSum.setAgentManager(agentManager);
    }
    
//--------------------------------------------------------------
// steps all particles across time interval tStep

    public void doStepInternal() {
        super.doStepInternal();
        //Compute forces on each atom
        forceSum.reset();
        potential.calculate(box, allAtoms, forceSum);
        pressureTensor.E(forceSum.getPressureTensor());

        //take step
        AtomSet leafList = box.getSpeciesMaster().getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            pressureTensor.E(forceSum.getPressureTensor());
            IVector v = a.getVelocity();
            workTensor.Ev1v2(v,v);
            workTensor.TE(((AtomTypeLeaf)a.getType()).getMass());
            pressureTensor.PE(workTensor);
            
            Agent agent = (Agent)agentManager.getAgent(a);
            IVector r = a.getPosition();
            work.E(r);
            r.PE(agent.rMrLast);
            agent.force.TE(((AtomTypeLeaf)a.getType()).rm()*t2);
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
    
    public void reset() throws ConfigurationOverlapException {
        super.reset();
        updateMrLast();
    }
    
    protected void updateMrLast() {
        AtomSet leafList = box.getSpeciesMaster().getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            Agent agent = (Agent)agentManager.getAgent(a);
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
        Agent agent = (Agent)agentManager.getAgent(atom);
        agent.rMrLast.Ea1Tv1(timeStep,atom.getVelocity());//06/13/03 removed minus sign before timeStep
    }
    
//--------------------------------------------------------------
    
    public Class getAgentClass() {
        return Agent.class;
    }

    public final Object makeAgent(IAtom a) {
        return new Agent(potential.getSpace());
    }
    
    public void releaseAgent(Object agent, IAtom atom) {}
            
	public final static class Agent implements IntegratorBox.Forcible {  //need public so to use with instanceof
        public IVector force;
        public IVector rMrLast;  //r - rLast

        public Agent(Space space) {
            force = space.makeVector();
            rMrLast = space.makeVector();
        }
        
        public IVector force() {return force;}
    }//end of Agent
    
}//end of IntegratorVerlet

