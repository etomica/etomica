package etomica;

/**
 * Evaluates the virial summed over all iterated atoms.
 *
 * @author David Kofke
 */
 
 /* History
  * 10/12/02 (DAK) new
  */
public class PotentialCalculationVirialSum extends PotentialCalculation 
											 implements PotentialCalculation.Summable {
											 	
	private Potential2.Soft p2Soft;
	
    protected double sum = 0.0;
        
    public PotentialCalculationVirialSum() {}
        
    public PotentialCalculation.Summable reset() {sum = 0.0; return this;}
    public double sum() {return sum;}

	public PotentialCalculation set(Potential2 p2) {
		if(!(p2 instanceof Potential2.Soft)) throw new RuntimeException("Error: PotentialCalculationVirialSum being used with potential that is not soft 2-body type");
		p2Soft = (Potential2.Soft)p2;
		return super.set(p2);
	}
    
    //pair
    public void actionPerformed(AtomPair pair) {
         sum += p2Soft.virial(pair);
    }//end of calculate
    
 }//end VirialSum
