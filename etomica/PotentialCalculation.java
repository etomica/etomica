package etomica;

/* History
 * 08/29/03 (DAK) made actionPerformed(AtomSet) abstract; added PotentialN
 */
public interface PotentialCalculation extends AtomsetActive {
 	
	public void doCalculation(AtomsetIterator iterator, Potential potential);
	
    public interface Summable {
        public double sum();
        public Summable reset();
    }
}//end of PotentialCalculation