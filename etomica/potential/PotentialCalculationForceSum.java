package etomica.potential;

import etomica.atom.Atom;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.iterator.AtomsetIterator;
import etomica.integrator.Integrator;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.
 */
public class PotentialCalculationForceSum extends PotentialCalculation {
        
    private final Vector f;
    protected Integrator.Forcible[] integratorAgents;
    
    public PotentialCalculationForceSum(Space space) {
        f = space.makeVector();
    }

    public void setAgents(Integrator.Forcible[] agents) {
        integratorAgents = agents;
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
					integratorAgents[((Atom)atoms).getGlobalIndex()].force().ME(f);
					break;
				case 2:
			        integratorAgents[((AtomPair)atoms).atom0.getGlobalIndex()].force().PE(f);
                    integratorAgents[((AtomPair)atoms).atom1.getGlobalIndex()].force().ME(f);
			 		break;
			 	default:
			 		throw new RuntimeException("Force calculation not implemented to treat (n>2)-body interactions");	
			}
		}
	}
}//end ForceSums
    
