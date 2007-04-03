package etomica.threaded;

import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.AtomsetIterator;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorPhase;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.phase.Phase;
import etomica.potential.Potential;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialCalculationForceSum;
import etomica.space.IVector;

public class PotentialCalculationForceSumThreaded extends PotentialCalculationForceSum implements PotentialCalculationThreaded, AgentSource{

	final protected PotentialCalculationForceSum[] pc;
	protected AtomAgentManager[] atomAgentManager;
    
	public PotentialCalculationForceSumThreaded(PotentialCalculationForceSum[] pc) {
		this.pc = pc;
	}

    public void reset(){
        super.reset();
        for (int i=0; i<pc.length; i++){
            pc[i].reset();
        }
    }
    
	public void setAgentManager(AtomAgentManager agentManager) {
        super.setAgentManager(agentManager);
        atomAgentManager = new AtomAgentManager[pc.length];
        
        for (int i=0; i<pc.length; i++){
            atomAgentManager[i] = new AtomAgentManager(this, agentManager.getPhase());
            pc[i].setAgentManager(atomAgentManager[i]);
            agentManager.getPhase();
		}
		
	}
	
	public void doCalculation(AtomsetIterator iterator, Potential potential) {
		throw new RuntimeException("This is not the correct 'doCalculation' to call.");
	}
	
	/* (non-Javadoc)
	 * @see etomica.threads.PotentialCalculationThreaded#getPotentialCalculations()
	 */
	public PotentialCalculation[] getPotentialCalculations(){
		return pc;
	}
	
	public void writeData(){
       
		Phase phase = integratorAgentManager.getPhase();
        AtomArrayList atomArrayList = phase.getSpeciesMaster().getLeafList();
      
        for(int j=0; j<atomArrayList.size(); j++){
            IVector force = ((IntegratorPhase.Forcible)integratorAgentManager.getAgent(atomArrayList.get(j))).force();
      
            for(int i=0; i<pc.length; i++){
                force.PE(((IntegratorPhase.Forcible)atomAgentManager[i].getAgent(atomArrayList.get(j))).force());
               
                
            }
        }
            
	}
    
    public Class getAgentClass() {
        return MyAgent.class;
    }

    public final Object makeAgent(Atom a) {
        return new MyAgent(integratorAgentManager.getPhase().getSpace());
    }
    
    public void releaseAgent(Object object, Atom atom){
        
    }
    
}
