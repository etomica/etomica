package etomica;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.
 */
public class PotentialCalculationForceSum implements PotentialCalculation {
        
    private final Space.Vector f;
    
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
			f.E(potentialSoft.gradient(atoms);
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
    
