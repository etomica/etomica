package etomica;

/**
 * @author David Kofke
 *
 * Parent interface for all iterators that return atoms, atom pairs, etc.
 */
public interface AtomSetIterator {

	public void all(AtomSet basis, IteratorDirective id, AtomSetAction action);
	
}
