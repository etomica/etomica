package etomica;

/**
 * Evaluates the virial summed over all iterated atoms.
 *
 * @author David Kofke
 */
 
 /* History
  * 10/12/02 (DAK) new
  * 08/29/03 (DAK) added actionPerformed(AtomSet) method because method made
  * abstract in PotentialCalculation
  */
public class PotentialCalculationVirialSum extends PotentialCalculation {
											 	
	/**
	 * Adds to the virial sum the energy values obtained from application of the given potential to the
	 * atoms produced by the given iterator.  Iterator is reset by method before beginning calculation.
	 */
	protected void doCalculation(AtomsetIterator iterator, Potential potential) {
		iterator.reset();
		while(iterator.hasNext()) {
			sum += ((Potential2.Soft)potential).virial(iterator.next());
		}
	}
	
	/**
	 * Sets the virial sum to zero, typically to begin a new virial-sum calculation.
	 * @return this instance, so the method can be called in-line as the instance is
	 * passed to the PotentialMaster.
	 */
	public PotentialCalculationVirialSum zeroSum() {
		sum = 0.0;
		return this;
	}

	/**
	 * Returns the current value of the energy sum.
	 */
	public double getSum() {return sum;}
	
	private double sum = 0.0;

 }//end VirialSum
