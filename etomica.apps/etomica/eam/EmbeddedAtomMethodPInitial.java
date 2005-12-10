package etomica.potential;
import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.phase.Phase;
import etomica.potential.EmbeddedAtomMethodP2.Wrapper;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * This class initializes the allatomAgents[EmbeddedAtomMethodP2.agentIndex] array.
 * All atoms' energy values and gradients are set to 0 and (0,0,0), respectively.
 * 
 * EmbeddedAtomMethodPInitial is a one-body potential, and is used along with 
 * EmbeddedAtomMethodP2 and EmbeddedAtomMethodPMany in the EAMMd3D simulation class.
 * 
 * This class was created by A. Schultz and K.R. Schadel July 2005.
 */

public final class EmbeddedAtomMethodPInitial extends Potential1 implements PotentialSoft {

	public EmbeddedAtomMethodPInitial(Space space, EmbeddedAtomMethodP2 eamP2) {
		super(space);
		gradient = space.makeVector();
		gradient.E(0);
        agentSource = eamP2;
    }
	
    public void setPhase(Phase phase) {
        super.setPhase(phase);
        agents = agentSource.getAgents(phase);
    }
    
	public double energy(AtomSet a) {
		agents[((Atom)a).getGlobalIndex()].x = 0;
		return 0;
	}

	public Vector gradient(AtomSet a) {
		agents[((Atom)a).getGlobalIndex()].x = 0;
        agents[((Atom)a).getGlobalIndex()].A.E(0);
		return gradient;
	}
	
	public double virial(AtomSet atoms) {
	    return 0.0;
    }
   
	Vector gradient;
    private final EmbeddedAtomMethodP2 agentSource;
    private Wrapper[] agents;
}
