package etomica.potential;
import etomica.atom.Atom;
import etomica.atom.AtomSet;
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

public final class EmbeddedAtomMethodPInitial extends Potential1 {

	public EmbeddedAtomMethodPInitial(Space space) {
		super(space);
		gradient = space.makeVector();
		gradient.E(0);
    }
	
	public double energy(AtomSet a) {
		((Wrapper)((Atom)a).allatomAgents[EmbeddedAtomMethodP2.agentIndex]).x = 0;
		return 0;
	}

	public Vector gradient(AtomSet a) {
		((Wrapper)((Atom)a).allatomAgents[EmbeddedAtomMethodP2.agentIndex]).x = 0;
		((Wrapper)((Atom)a).allatomAgents[EmbeddedAtomMethodP2.agentIndex]).A.E(0);
		return gradient;
	}
	
	public double virial(AtomSet atoms) {
	    return 0.0;
    }
   
	Vector gradient;
}

    

