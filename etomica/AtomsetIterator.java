package etomica;

/**
 * @author David Kofke
 *
 * Parent interface for all iterators that return atoms, atom pairs, etc.
 */
public interface AtomSetIterator {

	public void all(AtomSet basis, IteratorDirective id, AtomSetActive action);
	
	/**
	 * No-op instance of iterator that performs no action.
	 */
	public static final AtomSetIterator NULL = new AtomSetIterator() {
		public void all(AtomSet basis, IteratorDirective id, AtomSetActive action) {}
	};
}
