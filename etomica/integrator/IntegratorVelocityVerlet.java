package etomica.integrator;

import java.io.Serializable;

import etomica.EtomicaInfo;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.IAtom;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.exception.ConfigurationOverlapException;
import etomica.phase.Phase;
import etomica.potential.PotentialCalculationForcePressureSum;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.ICoordinateKinetic;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.util.IRandom;

public class IntegratorVelocityVerlet extends IntegratorMD implements AgentSource {

    private static final long serialVersionUID = 2L;
    protected PotentialCalculationForceSum forceSum;
    private final IteratorDirective allAtoms;
    protected final Tensor pressureTensor;
    protected final Tensor workTensor;
    
    protected AtomAgentManager agentManager;

    public IntegratorVelocityVerlet(Simulation sim) {
        this(sim.getPotentialMaster(),sim.getRandom(),
             sim.getDefaults().timeStep,sim.getDefaults().temperature);
    }
    
    public IntegratorVelocityVerlet(PotentialMaster potentialMaster, IRandom random,
            double timeStep, double temperature) {
        super(potentialMaster,random,timeStep,temperature);
        // if you're motivated to throw away information earlier, you can use 
        // PotentialCalculationForceSum instead.
        forceSum = new PotentialCalculationForcePressureSum(potentialMaster.getSpace());
        allAtoms = new IteratorDirective();
        // allAtoms is used only for the force calculation, which has no LRC
        // but we're also calculating the pressure tensor, which does have LRC.
        // things deal with this OK.
        allAtoms.setIncludeLrc(true);
        pressureTensor = potentialMaster.getSpace().makeTensor();
        workTensor = potentialMaster.getSpace().makeTensor();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Molecular dynamics using velocity Verlet integration algorithm");
        return info;
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

    // assumes one phase
    public void doStepInternal() {
        AtomArrayList leafList = phase.getSpeciesMaster().getLeafList();
        int nLeaf = leafList.size();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            AtomLeaf a = (AtomLeaf)leafList.get(iLeaf);
            MyAgent agent = (MyAgent)agentManager.getAgent(a);
            IVector r = a.getPosition();
            IVector v = ((ICoordinateKinetic)a).getVelocity();
            v.PEa1Tv1(0.5*timeStep*((AtomTypeLeaf)a.getType()).rm(),agent.force);  // p += f(old)*dt/2
            r.PEa1Tv1(timeStep,v);         // r += p*dt/m
        }

        forceSum.reset();
        //Compute forces on each atom
        potential.calculate(phase, allAtoms, forceSum);
        
        if(forceSum instanceof PotentialCalculationForcePressureSum){
            pressureTensor.E(((PotentialCalculationForcePressureSum)forceSum).getPressureTensor());
        }
        
        //Finish integration step
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            AtomLeaf a = (AtomLeaf)leafList.get(iLeaf);
//            System.out.println("force: "+((MyAgent)a.ia).force.toString());
            IVector velocity = ((ICoordinateKinetic)a).getVelocity();
            workTensor.Ev1v2(velocity,velocity);
            workTensor.TE(((AtomTypeLeaf)a.getType()).getMass());
            pressureTensor.PE(workTensor);
            velocity.PEa1Tv1(0.5*timeStep*((AtomTypeLeaf)a.getType()).rm(),((MyAgent)agentManager.getAgent(a)).force);  //p += f(new)*dt/2
        }
        
        pressureTensor.TE(1/phase.getBoundary().volume());

        if(isothermal) {
            doThermostat();
        }
        return;
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

        forceSum.reset();
        potential.calculate(phase, allAtoms, forceSum);
    }
              
//--------------------------------------------------------------
    
    public Class getAgentClass() {
        return MyAgent.class;
    }

    public final Object makeAgent(IAtom a) {
        return new MyAgent(potential.getSpace());
    }
    
    public void releaseAgent(Object agent, IAtom atom) {}
            
    public final static class MyAgent implements IntegratorPhase.Forcible, Serializable {  //need public so to use with instanceof
        private static final long serialVersionUID = 1L;
        public IVector force;

        public MyAgent(Space space) {
            force = space.makeVector();
        }
        
        public IVector force() {return force;}
    }
    
}
