package etomica.integrator;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.Atom.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.exception.ConfigurationOverlapException;
import etomica.phase.Phase;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.ICoordinateKinetic;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Constant NVT Molecular Dynamics Integrator-Constraint Method
 * 
 * Uses a modified version of the Leap Frog Algorithm. 
 * ke is calculated from the unconstrained velocities at time T.
 * A ratio of the setTemperature to the unconstrained temp (as solved from ke),
 * is used to calculate the new constrained velocities at T+Dt/2.  
 * The positions at T+Dt are solved for from the constrained velocities calculated at T+Dt/2.
 *  
 * @author Chris Iacovella
 * @author David Kofke
 */
public final class IntegratorConNVT extends IntegratorMD implements EtomicaElement, AgentSource {

    public final PotentialCalculationForceSum forceSum;
    private final IteratorDirective allAtoms = new IteratorDirective();
    private final Space space;
    Vector work, work1, work2, work3, work4;
    double halfTime, mass;

    protected Agent[] agents;
    protected AtomAgentManager agentManager;

    public IntegratorConNVT(Simulation sim) {
        this(sim.potentialMaster,sim.space,sim.getDefaults().timeStep,sim.getDefaults().temperature);
    }
    
    public IntegratorConNVT(PotentialMaster potentialMaster, Space space, 
            double timeStep, double temperature) {
        super(potentialMaster,timeStep,temperature);
        forceSum = new PotentialCalculationForceSum(space);
        this.space = space;
        work = space.makeVector();
        work1 = space.makeVector();
        work2 = space.makeVector();
        work3 = space.makeVector();
       	work4 = space.makeVector();
    }

	
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("NVT-MD Integrator using modified Leap-Frog algorithm");
        return info;
    }

	public boolean addPhase(Phase p) {
	    if(!super.addPhase(p)) return false;
        agentManager = new AtomAgentManager(this,p);
        return true;
    }
    
    public void removePhase(Phase oldPhase) {
        if (oldPhase == firstPhase) {
            // allow agentManager to de-register itself as a PhaseListener
            agentManager.setPhase(null);
            agentManager = null;
            agents = null;
        }
        super.removePhase(oldPhase);
    }
    
  	public final void setTimeStep(double t) {
    	super.setTimeStep(t);
    	halfTime = timeStep/2.0;
  	}
  	
  	private double Temper;
	public void setTemp(double temperature) {
		   Temper=temperature;
	   }
          
//--------------------------------------------------------------
// steps all particles across time interval tStep

    public void doStep() {
		double dim = space.D();  //get the dimension
		
        //Compute forces on each atom
        atomIterator.reset();
        while(atomIterator.hasNext()) {   //zero forces on all atoms
            agents[atomIterator.nextAtom().getGlobalIndex()].force.E(0.0);
        }
        potential.calculate(firstPhase, allAtoms, forceSum);
	
		//MoveA
		//Advance velocities from T-Dt/2 to T without constraint
		atomIterator.reset();	
		double Free=0.0;
		//degrees of freedom
		Free=((firstPhase.moleculeCount()-1)*dim); 
		
		double k=0.0;
        double chi;
		while(atomIterator.hasNext()) {
			Atom a = atomIterator.nextAtom();
			Vector v = ((ICoordinateKinetic)a.coord).velocity();
            
			work1.E(v); //work1 = v
			work2.E(agents[a.getGlobalIndex()].force);	//work2=F
			work1.PEa1Tv1(halfTime*((AtomTypeLeaf)a.type).rm(),work2); //work1= p/m + F*Dt2/m = v + F*Dt2/m
            
        	k+=work1.squared();
 
		}   
    	//calculate scaling factor chi
		
		chi= Math.sqrt(Temper*Free/(mass*k));
		
		//calculate constrained velocities at T+Dt/2
		atomIterator.reset();
		while(atomIterator.hasNext()) {
			Atom a = atomIterator.nextAtom();
			Agent agent = agents[a.getGlobalIndex()];
			Vector v = ((ICoordinateKinetic)a.coord).velocity();
		
			double scale = (2.0*chi-1.0); 
			work3.Ea1Tv1(scale,v); 
			work4.Ea1Tv1(chi*((AtomTypeLeaf)a.type).rm(),agent.force);
			work4.TE(timeStep);
			work3.PE(work4);
			v.E(work3);
				
		} 

		atomIterator.reset();
		while(atomIterator.hasNext()) {
			Atom a = atomIterator.nextAtom();
			Vector r = a.coord.position();
			Vector v = ((ICoordinateKinetic)a.coord).velocity();
            
			work.Ea1Tv1(timeStep,v);
			work.PE(r);
			r.E(work);
		}

        
    }//end of doStep
    

    public void reset() throws ConfigurationOverlapException {
        // reset might be called because atoms were added or removed
        // calling getAgents ensures we have an up-to-date array.
        agents = (Agent[])agentManager.getAgents();
        forceSum.setAgents(agents);
        super.reset();
    }

    public final Object makeAgent(Atom a) {
        return new Agent(space,a);
    }
            
	public final static class Agent implements Integrator.Forcible {  //need public so to use with instanceof
        public Atom atom;
        public Vector force;
    

        public Agent(Space space, Atom a) {
            atom = a;
            force = space.makeVector();
        }
        
        public Vector force() {return force;}
    }//end of Agent
    
}

