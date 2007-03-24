package etomica.integrator;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.phase.Phase;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.ICoordinateKinetic;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.util.IRandom;

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

    private static final long serialVersionUID = 1L;
    public final PotentialCalculationForceSum forceSum;
    private final IteratorDirective allAtoms;
    IVector work, work1, work2, work3, work4;
    double halfTime, mass;

    protected AtomAgentManager agentManager;

    public IntegratorConNVT(Simulation sim) {
        this(sim.getPotentialMaster(),sim.getRandom(),sim.getDefaults().timeStep,sim.getDefaults().temperature);
    }
    
    public IntegratorConNVT(PotentialMaster potentialMaster, IRandom random, 
            double timeStep, double temperature) {
        super(potentialMaster,random,timeStep,temperature);
        forceSum = new PotentialCalculationForceSum();
        allAtoms = new IteratorDirective();
        // allAtoms is used only for the force calculation, which has no LRC
        allAtoms.setIncludeLrc(false);
        Space space = potentialMaster.getSpace();
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

    public void setPhase(Phase p) {
        if (phase != null) {
            agentManager.dispose();
        }
        super.setPhase(p);
        agentManager = new AtomAgentManager(this,p);
        forceSum.setAgentManager(agentManager);
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
		double dim = phase.getSpace().D();  //get the dimension
		
        //Compute forces on each atom
        atomIterator.reset();
        while(atomIterator.hasNext()) {   //zero forces on all atoms
            ((Agent)agentManager.getAgent(atomIterator.nextAtom())).force.E(0.0);
        }
        potential.calculate(phase, allAtoms, forceSum);
	
		//MoveA
		//Advance velocities from T-Dt/2 to T without constraint
		atomIterator.reset();	
		double Free=0.0;
		//degrees of freedom
		Free=((phase.moleculeCount()-1)*dim); 
		
		double k=0.0;
        double chi;
		while(atomIterator.hasNext()) {
			AtomLeaf a = (AtomLeaf)atomIterator.nextAtom();
			IVector v = ((ICoordinateKinetic)a.getCoord()).getVelocity();
            
			work1.E(v); //work1 = v
			work2.E(((Agent)agentManager.getAgent(a)).force);	//work2=F
			work1.PEa1Tv1(halfTime*((AtomTypeLeaf)a.getType()).rm(),work2); //work1= p/m + F*Dt2/m = v + F*Dt2/m
            
        	k+=work1.squared();
 
		}   
    	//calculate scaling factor chi
		
		chi= Math.sqrt(Temper*Free/(mass*k));
		
		//calculate constrained velocities at T+Dt/2
		atomIterator.reset();
		while(atomIterator.hasNext()) {
			AtomLeaf a = (AtomLeaf)atomIterator.nextAtom();
			Agent agent = (Agent)agentManager.getAgent(a);
			IVector v = ((ICoordinateKinetic)a.getCoord()).getVelocity();
		
			double scale = (2.0*chi-1.0); 
			work3.Ea1Tv1(scale,v); 
			work4.Ea1Tv1(chi*((AtomTypeLeaf)a.getType()).rm(),agent.force);
			work4.TE(timeStep);
			work3.PE(work4);
			v.E(work3);
				
		} 

		atomIterator.reset();
		while(atomIterator.hasNext()) {
			AtomLeaf a = (AtomLeaf)atomIterator.nextAtom();
			IVector r = a.getCoord().getPosition();
			IVector v = ((ICoordinateKinetic)a.getCoord()).getVelocity();
            
			work.Ea1Tv1(timeStep,v);
			work.PE(r);
			r.E(work);
		}

        
    }//end of doStep
    

    public Class getAgentClass() {
        return Agent.class;
    }
    
    public final Object makeAgent(Atom a) {
        return new Agent(phase.getSpace());
    }
    
    public void releaseAgent(Object agent, Atom atom) {}
            
	public final static class Agent implements IntegratorPhase.Forcible {  //need public so to use with instanceof
        public IVector force;

        public Agent(Space space) {
            force = space.makeVector();
        }
        
        public IVector force() {return force;}
    }
    
}
