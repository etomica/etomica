package etomica.action;

import etomica.Action;
import etomica.Atom;

/**
 * Interface for a class that can perform an action on an atom set.
 */

/* History
 * 
 * 01/25/03 (DAK) new
 */
public interface AtomsetAction extends Action {
    
	public void actionPerformed(Atom[] atoms);
    public void setAtoms(Atom[] atom);
    public Atom[] getAtoms();
    
    public static final AtomsetAction NULL = new AtomsetActionAdapter() {
        public final void actionPerformed(Atom[] a) {}
    };

}
