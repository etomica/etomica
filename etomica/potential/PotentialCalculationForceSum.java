package etomica.potential;

import etomica.atom.Atom;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.iterator.AtomsetIterator;
import etomica.integrator.IntegratorPhase;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.
 */
public class PotentialCalculationForceSum extends PotentialCalculation {
        
    private Vector[] f;
    protected IntegratorPhase.Forcible[] integratorAgents;
    
    public PotentialCalculationForceSum(Space space) {
//        f = space.makeVector();
    }

    public void setAgents(IntegratorPhase.Forcible[] agents) {
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
			f = potentialSoft.gradient(atoms);
	//TODO update gradient method to return array of gradient vectors, one for each atom
	//TODO make a consistent definition so that atoms[0] is always ME(f)
			switch(nBody) {
				case 1: 
					integratorAgents[((Atom)atoms).getGlobalIndex()].force().ME(f[0]);
					break;
				case 2:
			        integratorAgents[((AtomPair)atoms).atom0.getGlobalIndex()].force().ME(f[0]);
                    integratorAgents[((AtomPair)atoms).atom1.getGlobalIndex()].force().ME(f[1]);
			 		break;
                default:
                    //XXX atoms.count might not equal f.length.  The potential is responsible for
                    // returning an array that includes force vectors for atoms it determined to be
                    // non-zero.  It's allowed to omit the vectors for Atoms that are 0.
                    for (int i=0; i<f.length; i++) {
                        integratorAgents[atoms.getAtom(i).getGlobalIndex()].force().ME(f[i]);
                    }
			}
		}
	}
}//end ForceSums
    
