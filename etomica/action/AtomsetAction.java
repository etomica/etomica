package etomica.action;

import etomica.Atom;

/**
 * @author David Kofke
 *
 * Interface for a class that can perform an action on an atom set.
 */

/* History
 * 
 * 01/25/03 (DAK) new
 */
public interface AtomsetAction {
    
	public void actionPerformed(Atom[] atoms);
    public void setAtoms(Atom[] atom);
    public Atom[] getAtoms();
    
    public static final AtomsetAction NULL = new AtomsetActionAdapter() {
        public final void actionPerformed(Atom[] a) {}
    };

}
