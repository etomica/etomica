package etomica;

/* History
 * 08/29/03 (DAK) made actionPerformed(AtomSet) abstract; added PotentialN
 */


public interface PotentialCalculation {
 	
	public void doCalculation(AtomsetIterator iterator, Potential potential);
	
}//end of PotentialCalculation