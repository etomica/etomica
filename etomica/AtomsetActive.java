package etomica;

/**
 * @author David Kofke
 *
 * Interface for a class that can perform an action on an atom set.
 */

/* History
 * 
 * 01/25/03 (DAK) new
 */
public interface AtomsetActive {
	public void actionPerformed(Atom[] atoms);
}
