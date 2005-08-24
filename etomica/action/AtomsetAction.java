package etomica.action;

import etomica.AtomSet;

/**
 * Interface for a class that can perform an action on an atom set.
 */

/* History
 * 
 * 01/25/03 (DAK) new
 */
public interface AtomsetAction extends Action {
    
	public void actionPerformed(AtomSet atoms);
    public void setAtoms(AtomSet atom);
    public AtomSet getAtoms();
    
    public static final AtomsetAction NULL = new AtomsetActionAdapter() {
        public final void actionPerformed(AtomSet a) {}
    };

}
