package etomica.integrator;

import java.io.Serializable;

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

/* History of changes
 * 08/29/02 (DAK) changed Andersen thermostat to velocity-scaling thermostat
 * 01/10/03 (DAK) reintroduced Andersen thermostat with flag to allow
 * selection of it or velocity-rescaling as the thermostat mechanism
 * */

public final class IntegratorVelocityVerlet extends IntegratorMD implements EtomicaElement, AgentSource {

    public final PotentialCalculationForceSum forceSum;
    private final Space space;
    private final IteratorDirective allAtoms;
    
    protected AtomAgentManager agentManager;

    public IntegratorVelocityVerlet(Simulation sim) {
        this(sim.potentialMaster,sim.space,sim.getDefaults().timeStep,sim.getDefaults().temperature);
    }
    
    public IntegratorVelocityVerlet(PotentialMaster potentialMaster, Space space,
            double timeStep, double temperature) {
        super(potentialMaster,timeStep,temperature);
        this.space = space;
        forceSum = new PotentialCalculationForceSum();
        allAtoms = new IteratorDirective();
        // allAtoms is used only for the force calculation, which has no LRC
        allAtoms.setIncludeLrc(false);
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
    }
    
//--------------------------------------------------------------
// steps all particles across time interval tStep

    // assumes one phase
    public void doStep() {
        atomIterator.setPhase(phase);
        atomIterator.reset();              //reset iterator of atoms
        while(atomIterator.hasNext()) {    //loop over all atoms
            AtomLeaf a = (AtomLeaf)atomIterator.nextAtom();  //  advancing positions full step
            MyAgent agent = (MyAgent)agentManager.getAgent(a);     //  and momenta half step
            Vector r = a.coord.position();
            Vector v = ((ICoordinateKinetic)a.coord).velocity();
            v.PEa1Tv1(0.5*timeStep*((AtomTypeLeaf)a.type).rm(),agent.force);  // p += f(old)*dt/2
            r.PEa1Tv1(timeStep,v);         // r += p*dt/m
            agent.force.E(0.0);
        }
                
        //Compute forces on each atom
        potential.calculate(phase, allAtoms, forceSum);
        
        //Finish integration step
        atomIterator.reset();
        while(atomIterator.hasNext()) {     //loop over atoms again
            AtomLeaf a = (AtomLeaf)atomIterator.nextAtom();   //  finishing the momentum step
//            System.out.println("force: "+((MyAgent)a.ia).force.toString());
            ((ICoordinateKinetic)a.coord).velocity().PEa1Tv1(0.5*timeStep*((AtomTypeLeaf)a.type).rm(),((MyAgent)agentManager.getAgent(a)).force);  //p += f(new)*dt/2
        }
        if(isothermal) {
            doThermostat();
        }
        return;
    }
    
    
//--------------------------------------------------------------

    public void reset() throws ConfigurationOverlapException{
        if(!initialized) return;
        // reset might be called because atoms were added or removed
        // calling getAgents ensures we have an up-to-date array.
        forceSum.setAgentManager(agentManager);
        
        atomIterator.setPhase(phase);
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
            ((MyAgent)agentManager.getAgent(a)).force.E(0.0);
        }
        potential.calculate(phase, allAtoms, forceSum);//assumes only one phase
        super.reset();
    }
              
//--------------------------------------------------------------
    
    public Class getAgentClass() {
        return MyAgent.class;
    }

    public final Object makeAgent(Atom a) {
        return new MyAgent(space,a);
    }
    
    public void releaseAgent(Object agent, Atom atom) {}
            
    public final static class MyAgent implements IntegratorPhase.Forcible, Serializable {  //need public so to use with instanceof
        public Atom atom;
        public Vector force;

        public MyAgent(Space space, Atom a) {
            atom = a;
            force = space.makeVector();
        }
        
        public Vector force() {return force;}
    }//end of MyAgent
    
}//end of IntegratorVelocityVerlet

