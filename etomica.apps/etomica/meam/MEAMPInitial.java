package etomica.meam;
import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomSet;
import etomica.meam.MEAMP2.Wrapper;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Vector3D;

/**
 * This class initializes the values of every element in the arrays "sums" and
 * "gradientSums" to be zero or the zero vector ("nada"), respectively.  
 * 
 * This class provides no direct contribution to the values of either the 
 * potential or the gradient of the potential, although it does function as 
 * a one-body potential.  Note that the energy() method and gradient() method 
 * return zero and the zero vector, respectively. 
 * 
 * This class was created by A. Schultz and K.R. Schadel July 2005 as part of
 * a pseudo embedded-atom method potential.  In February 2006 it was adapted 
 * to be a part of a modified embedded-atom method potential.
 */

public final class MEAMPInitial extends Potential1 implements PotentialSoft {

	public MEAMPInitial(Space space, PhaseAgentManager phaseAgentManager) {
		super(space);
		this.phaseAgentManager = phaseAgentManager;
    }
	
	public void setPhase(Phase phase) {
		AtomAgentManager[] agentManager = (AtomAgentManager[])phaseAgentManager.getAgents();
        agents = (Wrapper[])agentManager[phase.getIndex()].getAgents();
	}
    
	public double energy(AtomSet a) {
		Wrapper agent = agents[((Atom)a).getGlobalIndex()];
		for(int i=0; i < 27; i++) {
			agent.sums[i] = 0;
		}
		//MEAMPInitial should contribute nothing to the potential directly.
		return 0;
	}
	
	public Vector[] gradient(AtomSet a) {
		Wrapper agent = agents[((Atom)a).getGlobalIndex()];
		for(int i=0; i < 27; i++) {
			agent.sums[i] = 0; //energy(a) may not be called
			agent.gradientSums[i].E(0);
		}
		//MEAMPInitial should contribute nothing to the gradient directly.
		return nada; 
	}
	
	public double virial(AtomSet atoms) {
	    return 0.0;
    }
   
	Vector gradient;
    private Wrapper[] agents;
    private PhaseAgentManager phaseAgentManager;
	private final Vector3D[] nada = new Vector3D[0];
}
