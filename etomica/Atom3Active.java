package etomica;

/**
 * @author David Kofke
 *
 * Interface for a class that can perform an action on an atom triplet.
 */

/* History
 * 
 * 01/25/03 (DAK) new
 */
public interface Atom3Active extends AtomsetActive {
	public void actionPerformed(Atom3 atom3);
}
