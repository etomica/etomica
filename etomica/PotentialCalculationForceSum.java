package etomica;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.
 */
public final class PotentialCalculationForceSum implements Potential1Calculation, Potential2Calculation {
        
    private final Space.Vector f;
    public PotentialCalculationForceSum(Space space) {
            f = space.makeVector();
    }
        
    //atom
    public void calculate(AtomIterator iterator, Potential1 potential) {
        Potential1Soft potentialSoft = (Potential1Soft)potential;
        while(iterator.hasNext()) {
            Atom atom = iterator.next();
            f.E(potentialSoft.gradient(atom));
            ((Integrator.Agent.Forcible)atom.ia).force().ME(f);
        }//end while
    }//end of calculate

    //pair
    public void calculate(AtomPairIterator iterator, Potential2 potential) {
        Potential2Soft potentialSoft = (Potential2Soft)potential;
        while(iterator.hasNext()) {
            AtomPair pair = iterator.next();
            f.E(potentialSoft.gradient(pair));
            ((Integrator.Agent.Forcible)pair.atom1().ia).force().PE(f);
            ((Integrator.Agent.Forcible)pair.atom2().ia).force().ME(f);
        }//end while
    }//end of calculate
}//end ForceSums
    
