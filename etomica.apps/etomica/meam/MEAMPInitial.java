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
 * This class initializes the allatomAgents[EmbeddedAtomMethodP2.agentIndex] array.
 * All atoms' energy values and gradients are set to 0 and (0,0,0), respectively.
 * 
 * EmbeddedAtomMethodPInitial is a one-body potential, and is used along with 
 * EmbeddedAtomMethodP2 and EmbeddedAtomMethodPMany in the EAMMd3D simulation class.
 * 
 * This class was created by A. Schultz and K.R. Schadel July 2005.
 */

public final class MEAMPInitial extends Potential1 implements PotentialSoft {

	public MEAMPInitial(Space space, PhaseAgentManager phaseAgentManager) {
		super(space);
		nada = (Vector3D)space.makeVector();
		nada.E(0);
		this.phaseAgentManager = phaseAgentManager;
    }
	
	public void setPhase(Phase phase) {
        agentManager = (AtomAgentManager)phaseAgentManager.getAgents()[phase.getIndex()];
    }
    
	public double energy(AtomSet a) {
		Wrapper agent = agents[((Atom)a).getGlobalIndex()];
		for(int i=0; i < 27; i++) {
			agent.sums[i] = 0;
		}
		return 0;
	}
	
	public Vector gradient(AtomSet a) {
		Wrapper agent = agents[((Atom)a).getGlobalIndex()];
		for(int i=0; i < 27; i++) {
			agent.gradientSums[i].E(0);
		}
		return nada; 
	}
	
	public double virial(AtomSet atoms) {
	    return 0.0;
    }
   
	Vector gradient;
    private Wrapper[] agents;
    private PhaseAgentManager phaseAgentManager;
	private AtomAgentManager agentManager;
	private final Vector3D nada;
}
