package etomica.potential;

import etomica.Space;
import etomica.atom.Atom;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.iterator.AtomsetIterator;
import etomica.integrator.Integrator;
import etomica.space.Vector;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.
 * 
 * Intended for use in the Molecular-Dynamics simulation employing the Embedded-Atom 
 * Method potential, EAMMd3D.java.
 * 
 * Created by A. Schultz and K.R. Schadel July 2005.
 */
public class PotentialCalculationForceSumEAM extends PotentialCalculation {
        
    private final Vector f;
    
    public PotentialCalculationForceSumEAM(Space space) {
            f = space.makeVector();
    }
        
	/**
	 * Adds forces due to given potential acting on the atoms produced by the iterator.
	 * Implemented for only 1- and 2-body potentials.
	 */
	public void doCalculation(AtomsetIterator iterator, Potential potential) {
		PotentialSoft potentialSoft = (PotentialSoft)potential;
		int nBody = potential.nBody();
		iterator.reset();
		while(iterator.hasNext()) {
			AtomSet atoms = iterator.next();
			f.E(potentialSoft.gradient(atoms));
	//TODO update gradient method to return array of gradient vectors, one for each atom
	//TODO make a consistent definition so that atoms[0] is always ME(f)
			switch(nBody) {
				case 1: 
					((Integrator.Forcible)((Atom)atoms).ia).force().ME(f);
					break;
				case 2: 
			        ((Integrator.Forcible)((AtomPair)atoms).atom0.ia).force().PE(f);
			        ((Integrator.Forcible)((AtomPair)atoms).atom1.ia).force().ME(f);
			 		break;
			 	default:
			 		throw new RuntimeException("Force calculation not implemented to treat (n>2)-body interactions");	
			}
		}
	}
}//end ForceSums
    
