package etomica;

/**
 * @author David Kofke
 *
 * Interface for a class that can perform an action on an atom.
 */

/* History
 * 
 * 01/25/03 (DAK) new
 */
public interface AtomActive extends AtomSetActive {
	public void actionPerformed(Atom atom);
}
