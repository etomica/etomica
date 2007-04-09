package etomica.integrator;

import etomica.EtomicaInfo;
import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.exception.ConfigurationOverlapException;
import etomica.phase.Phase;
import etomica.potential.PotentialCalculationForcePressureSum;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.ICoordinateKinetic;
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

    protected AtomAgentManager agentManager;

    public IntegratorVerlet(Simulation sim) {
        this(sim.getPotentialMaster(),sim.getRandom(),sim.getDefaults().timeStep,
                sim.getDefaults().temperature);
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
          
    public void setPhase(Phase p) {
        if (phase != null) {
            // allow agentManager to de-register itself as a PhaseListener
            agentManager.dispose();
        }
        super.setPhase(p);
        agentManager = new AtomAgentManager(this,p);
        forceSum.setAgentManager(agentManager);
    }
    
//--------------------------------------------------------------
// steps all particles across time interval tStep

    public void doStepInternal() {

        //Compute forces on each atom       
        forceSum.reset();
        potential.calculate(phase, allAtoms, forceSum);
        pressureTensor.E(forceSum.getPressureTensor());

        //take step
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            AtomLeaf a = (AtomLeaf)atomIterator.nextAtom();
            pressureTensor.E(forceSum.getPressureTensor());
            IVector v = ((ICoordinateKinetic)a.getCoord()).getVelocity();
            workTensor.Ev1v2(v,v);
            workTensor.TE(((AtomTypeLeaf)a.getType()).getMass());
            pressureTensor.PE(workTensor);
            
            Agent agent = (Agent)agentManager.getAgent(a);
            IVector r = a.getCoord().getPosition();
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
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            AtomLeaf a = (AtomLeaf)atomIterator.nextAtom();
            Agent agent = (Agent)agentManager.getAgent(a);
            agent.rMrLast.Ea1Tv1(timeStep,((ICoordinateKinetic)a.getCoord()).getVelocity());//06/13/03 removed minus sign before timeStep
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
    protected void randomizeMomentum(AtomLeaf atom) {
        super.randomizeMomentum(atom);
        Agent agent = (Agent)agentManager.getAgent(atom);
        agent.rMrLast.Ea1Tv1(timeStep,((ICoordinateKinetic)atom).getVelocity());//06/13/03 removed minus sign before timeStep
    }
    
//--------------------------------------------------------------
    
    public Class getAgentClass() {
        return Agent.class;
    }

    public final Object makeAgent(Atom a) {
        return new Agent(potential.getSpace());
    }
    
    public void releaseAgent(Object agent, Atom atom) {}
            
	public final static class Agent implements IntegratorPhase.Forcible {  //need public so to use with instanceof
        public IVector force;
        public IVector rMrLast;  //r - rLast

        public Agent(Space space) {
            force = space.makeVector();
            rMrLast = space.makeVector();
        }
        
        public IVector force() {return force;}
    }//end of Agent
    
}//end of IntegratorVerlet

