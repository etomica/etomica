package etomica;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.
 */
public class PotentialCalculationForceSum extends PotentialCalculation {
        
    private final Space.Vector f;
    private Potential2.Soft p2Soft;
    private Potential1.Soft p1Soft;
    
    public PotentialCalculationForceSum(Space space) {
            f = space.makeVector();
    }
        
	public PotentialCalculation set(Potential1 p1) {
		if(!(p1 instanceof Potential1.Soft)) throw new RuntimeException("Error: PotentialCalculationForceSum being used with potential that is not soft 2-body type");
		p1Soft = (Potential1.Soft)p1;
		return super.set(p1);
	}
	public PotentialCalculation set(Potential2 p2) {
		if(!(p2 instanceof Potential2.Soft)) throw new RuntimeException("Error: PotentialCalculationForceSum being used with potential that is not soft 2-body type");
		p2Soft = (Potential2.Soft)p2;
		return super.set(p2);
	}
    //atom
    public void actionPerformed(Atom atom) {
        f.E(p1Soft.gradient(atom));
        ((Integrator.Agent.Forcible)atom.ia).force().ME(f);
    }//end of calculate

    //pair
    public void actionPerformed(AtomPair pair) {
        f.E(p2Soft.gradient(pair));
        ((Integrator.Agent.Forcible)pair.atom1().ia).force().PE(f);
        ((Integrator.Agent.Forcible)pair.atom2().ia).force().ME(f);
    }//end of calculate
    
    public void actionPerformed(Atom3 atom3) {
    	throw new etomica.exception.MethodNotImplementedException();
    }
}//end ForceSums
    
