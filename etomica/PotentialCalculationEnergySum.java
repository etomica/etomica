package etomica;

/**
 * Evaluates the energy summed over all iterated atoms.
 *
 * @author David Kofke
 */
public class PotentialCalculationEnergySum extends PotentialCalculation
											implements PotentialCalculation.Summable {
    protected double sum = 0.0;
        
    public PotentialCalculation.Summable reset() {sum = 0.0; return this;}
    public double sum() {return sum;}
    
 	public void actionPerformed(Phase phase) {
 		sum += potential0.energy(phase);
 	}
 	   
    public void actionPerformed(Atom atom) {
        sum += potential1.energy(atom);
    }
    
    public void actionPerformed(AtomPair pair) {
//    	System.out.println(pair.toString());
        sum += potential2.energy(pair);
    }
    
    public void actionPerformed(Atom3 atom3) {
    	sum += potential3.energy(atom3);
    }
        
}//end EnergySum
