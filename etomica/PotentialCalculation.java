package etomica;

/* History
 * 08/29/03 (DAK) made actionPerformed(AtomSet) abstract; added PotentialN
 */
public abstract class PotentialCalculation implements AtomsetActive {
 	
	public abstract void doCalculation(AtomsetIterator iterator, Potential potential);
	
    public interface Summable {
        public double sum();
        public Summable reset();
    }
}//end of PotentialCalculation