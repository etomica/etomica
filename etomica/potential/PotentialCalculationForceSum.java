package etomica.potential;

import etomica.Atom;
import etomica.AtomsetIterator;
import etomica.Integrator;
import etomica.Potential;
import etomica.Space;
import etomica.Integrator.Forcible;
import etomica.space.Vector;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.
 */
public class PotentialCalculationForceSum extends PotentialCalculation {
        
    private final Vector f;
    
    public PotentialCalculationForceSum(Space space) {
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
			Atom[] atoms = iterator.next();
			f.E(potentialSoft.gradient(atoms));
	//TODO update gradient method to return array of gradient vectors, one for each atom
	//TODO make a consistent definition so that atoms[0] is always ME(f)
			switch(nBody) {
				case 1: 
					((Integrator.Forcible)atoms[0].ia).force().ME(f);
					break;
				case 2: 
			        ((Integrator.Forcible)atoms[0].ia).force().PE(f);
			        ((Integrator.Forcible)atoms[1].ia).force().ME(f);
			 		break;
			 	default:
			 		throw new RuntimeException("Force calculation not implemented to treat (n>2)-body interactions");	
			}
		}
	}
}//end ForceSums
    
