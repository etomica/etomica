package etomica;

/**
 * @author kofke
 *
 * Interface for a set of atoms, such as Atom, AtomPair, etc.
 */
public interface AtomSet {

	public boolean contains(Atom a);
	
	public int atomCount();
	
}
